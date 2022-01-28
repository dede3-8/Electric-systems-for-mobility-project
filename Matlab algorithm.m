close all
clear all
clc

g   = 9.80665;


M=1610;
mu=0.012;
beta=0.05;
cx=0.32;
v0=[0
1.831828571
3.66320706
5.493235745
7.321016891
9.164342232
11.00359168
12.80049434
14.59165982
16.37624442
18.13473069
19.88507707
21.62649747
23.37690752
25.11675164
26.84527912
28.56175945
30.228112
31.88132466
33.52076495
].';


n=1;
for i=0:1:26      %this is the speed profile. It is computed so that it stays for nearly 10 minutes at 120km/h, it then decelerates to 100km/h in 20 seconds, stays at 100km/h for 10 minutes, and accelerates to 120km/h in 20 seconds. It does that 26 times.

  for k=20:1:600
    v0(k+1200*i)=33.52076495;
  endfor
  n=1;
  for k=601:1:620
    v0(k+1200*i)=33.52076495-n/3.6;
    n=n+1;
  endfor
  for k=620:1:1200
    v0(k+1200*i)=27.77777777777778;
  endfor 
n=1;
   for k=1201:1:1220
    v0(k+1200*i)=27.77777777777778+n/3.6;
    n=n+1;
  endfor
endfor


a=[1.831828571
1.831378489
1.830028684
1.827781146
1.843325341
1.839249448
1.796902661
1.791165482
1.784584598
1.758486271
1.750346381
1.741420394
1.750410051
1.739844126
1.728527482
1.71648032
1.666352551
1.653212661
1.639440292
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0].';
for i=0:1:26            %These is the acceleration profile in order to achieve that speed profile.
  for k=1:1:20
    a(k+(1+i)*1200)=1/3.6;
  endfor
  for k=20:1:600
    a(k+1200*i)=0;
  endfor
  for k=601:1:620
    a(k+1200*i)=-1/3.6;
  endfor
  for k=620:1:1200
    a(k+1200*i)=0;
  endfor
endfor

inc=[0
].';
for k=1:1:1450                %This is the slope profile from Spain, in which we divided the path into 10 different parts in order to linearize it. With those 10 points we calculated how many seconds it would take to cross that part and inserted them in here. In the case of Canada its just all values of 0 because Canada is mostly flat. 
  inc(k)=1.39503;
endfor
for k=1451:1:2573
  inc(k)=-0.49271137;
endfor
for k=2573:1:4354
  inc(k)=0.95036765;
endfor
for k=4354:1:12930
  inc(k)=-0.078244274;
endfor
for k=12930:1:15091
  inc(k)=0.61067;
endfor
for k=15091:1:16957
  inc(k)=-1.5807;
endfor
for k=16957:1:21737
  inc(k)=-0.053424657;
endfor
for k=21737:1:24912
  inc(k)=0.553608247;
endfor
for k=24912:1:27564
  inc(k)=-0.9061728;
endfor
for k=27564:1:31688
  inc(k)=0;
endfor

ti=[0:1:32419];    %This is the complete time of the 1000km divided into seconds

tf=[1:1:32420];    %This is also the complete time of the 1000km divided into seconds, but +1 in order to divide the path into small sections of 1 seconds

Eaerototal=0;       %These lines are just to initialize the energy consumptions in 0 to avoid errors.
Erollingtotal=0;
Einertiatotal=0;
Einertiatotalpos=0;
Einertiatotalneg=0;
Eaerototalpos=0;
Eaerototalneg=0;
Erollingtotalpos=0;
Erollingtotalneg=0;
Einctotal=0;
Pmax=0;
Etot=0;
Etotal(1)=0;
for ii=1:1:1914    %This is the main loop of the code. It goes from 1 to 31688, and each ii is a second in the real trip. Knowing this, instead of doing the whole trip in one loop, what we did was initialize the data for the whole trip, but then use the loops' starting value and finishing value in order to simulate the stops and know the energy consumption at each stop. So after each stop what we did was change the initial and end values of the loop to the ones of the next part of the trip until the next stop. 


    Paero=@(t)0.5*1.225*0.208*2.39*((v0(ii)+a(ii).*t).^3);        %Here we are calculating the power needed by each of the four forces as a function of time, and then what we do is integrate that function in order to calculate the energy consumption. But instead of doing this for the whole trip at a time, we have done it second by second using the ti and tf, using ti as the second before and tf as the second at the moment.

    q(ii)=integral(Paero,0,1);
    Eaerototal=Eaerototal+q(ii)/0.92;


    Prolling=@(t) M*g*mu.*(v0(ii)+a(ii).*t);
    j(ii)=integral(Prolling,0,1);
    Erollingtotal=Erollingtotal+j(ii)/0.92;
    
    Pinertia=@(t) M*(1+beta)*a(ii).*(v0(ii)+a(ii).*t);
    s(ii)=integral(Pinertia,0,1);
    Einertiatotal=Einertiatotal+s(ii)/0.92;
    
    
    Pinc=@(t) inc(ii)/100*M*g.*(v0(ii)+a(ii).*t);
    in(ii)=integral(Pinc,0,1);
    Einctotal=Einctotal+in(ii)/0.92;
    if s(ii)>=0                                 %In this part we calculate if the inertial power is either negative or positive, and according to that we use one efficiency or the other.
    Einertiatotalpos=Einertiatotalpos+s(ii);
    Erollingtotalpos=Erollingtotalpos+j(ii);
    Eaerototalpos=Eaerototalpos+q(ii);
    Etotal(ii)=Etot+(q(ii)+s(ii)+j(ii)+in(ii))/0.92;
    Etot=Etotal(ii);
    else
    Einertiatotalneg=Einertiatotalneg+s(ii);
    Erollingtotalneg=Erollingtotalneg+j(ii);
    Eaerototalneg=Eaerototalneg+q(ii);
    Etotal(ii)=Etotal(ii-1)+(q(ii)+s(ii)+j(ii)+in(ii))*0.92*0.92;
    end
    Pmax2=(Paero(1)+Pinertia(1)+Prolling(1)+Pinc(1))/0.92;        %This is the calculation of the maximum power for the path
    if Pmax2>=Pmax
        Pmax=Pmax2;
    end
    
    Finertia(ii)=M*(1+beta)*a(ii);          %Here we calculate the forces in each second. 
    Frolling(ii)=M*g*mu;
    Faero(ii)=0.5*1.225*0.208*2.39*((v0(ii)+a(ii)*1).^2);
    
if (v0(ii)&&a(ii))~=0           %Here we set the rolling resistance to 0 if the car is stopped, and to that value if it is running.
        Frolling(ii)=M*g*mu;
    else
        Frolling(ii)=0;
    end
    Fresist(ii)=Faero(ii)+Finertia(ii)+Frolling(ii);
    
    
end

%for jj=1:10:1914         %this part is in comments, because we don't need it for our results, but we have used it to plot forces, power and energy a lot of times in order to get a wider scope of the situation of the car. Now we are not using it in order to make the code run faster.
   % Eto=Etotal(jj);
   % te=tf(jj);
  %  figure(2);
   % plot(te,Eto,'g*')
   % xlabel('t(s)');
    %ylabel('Energy consumed (J)');
    %hold on
%end
Einertiatotalpos              %Here we display the values of each energy consumption created by the different forces. The two last values are the maximum power and the total energy consumption. In the case of the last value which is the total nergy consumption we have put the end between the parenthesis in order for it to show the last value, regardless of in which second of the path it is.
Einertiatotalneg
Einertiatotal
Eaerototalpos
Eaerototalneg
Eaerototal
Erollingtotalpos
Erollingtotalneg
Erollingtotal
Einctotal
Pmax
Etotal(end)