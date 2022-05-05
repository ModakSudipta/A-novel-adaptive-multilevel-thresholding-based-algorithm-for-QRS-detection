close all;
clear all;clc;
 
%% Loading the signal and setting Parameters
tic;
x=rdsamp('104.dat');
flag=1;
x1=x(:,1);
x2=x(:,2); 
Sa=x2;
%Sa=[x1(1:18000);x1(20000:87500); x1(101000:555000); x1(590000:650000)]; 
%% Segmentation Initialization
Nnaught=length(Sa);
if Nnaught<=50000
    segcoef=2;
elseif Nnaught<=150000
    segcoef=6;
elseif Nnaught<=500000
    segcoef=20;
elseif Nnaught<=650000
    segcoef=26;
elseif Nnaught<=1000000
    segcoef=40;
 elseif Nnaught<=1500000
    segcoef=60;
 elseif Nnaught<=1800000
    segcoef=72;
else
    segcoef=74;
end
N=floor(Nnaught/segcoef);
%% Sampling Frequency Initialization
if flag==1
    Fs=360;
    dec=150;
    coef=0.7; 
    minsep=0.32;
    prep=1.5;
    obtimum=60;
elseif flag==2
Fs = 250;
dec=100;
coef=1.25;
prep=1.5;
minsep=0.32;
obtimum=72;
elseif flag==3
    Fs = 250;
dec=110;
coef=0.7;
prep=1.5;
minsep=0.32;
obtimum=72;
else
    Fs=360;
    dec=100;
    coef=0.7;
    Sa=x1(1:500000);
end
 
coef1=0.4;
aperture=0.48;
 
n=length(Sa)/N;
j=0;
 
%% Peaks Detection Initialization
[Sc baseline]=removebaseline(Sa,Fs,0.1);
Sc=abs(Sc);
[Sd baseline]=removebaseline(Sa,Fs,0.05);
Sd=abs(Sd);
mu=[];
sigma=[];
pola1=[];
pola2=[];


 [pid lid] = findpeaks(Sc,'MINPEAKDISTANCE',round(0.2*Fs));
 value1=mean(pid);
 [pid1 lid1] = findpeaks(Sc,'MINPEAKDISTANCE',round(0.2*Fs),'MINPEAKHEIGHT',coef1*value1);
 
 pola1=[pola1 lid1];
 
[pid lid] = findpeaks(Sd,'MINPEAKDISTANCE',round(0.2*Fs));
  value2=mean(pid);
  [pid2 lid2] = findpeaks(Sd,'MINPEAKDISTANCE',round(0.2*Fs),'MINPEAKHEIGHT',coef1*value2);
  
decision=length(pid1)*100/length(pid2);
if decision>dec
    multi=0.05;
else
    multi=0.1;
end
if flag==1
gror=rdsamp('231.dat');
gror1=gror(:,1);
gror4=rdsamp('111.dat');
gror5=gror4(:,1);
% if Sa==gror1 | Sa==gror5
%     multi=0.05;
% % elseif Sa==-1*gior1
% %     coef=1;
% %     multi=0.5;
% else
% end
% else
 end
 


%% Preprocessing using two median filters and moving mean algorithm
[SS baseline]=removebaseline(Sa,Fs,multi);
 
for l=1:length(SS)
    if l>10 && l<length(SS)-10
        s1(l)=sum(SS(l-10:l+9));
        avg1=s1(l)/20;
        Sb(l)=avg1;
    else
        Sb(l)=SS(l);
    end
end

Sc=Sb; 
Sb=abs(Sb);
 
%% False Peak elimination due to high frequency noise
n=floor(length(SS)/N);
j=0;
 
locations=[];
peaks=[];
finalpeaks=[];
finallocations=[];
refv=[];
counterloc=[];
refstd=[];
for i=1:n
 
        k=i*N;
        y=Sb(j+1:k);
        z=Sa(j+1:k);
        [pks1,locs1] = findpeaks(y,'MINPEAKDISTANCE',round(0.2*Fs));
 
         value=mean(pks1);
         if flag==2
         if i>1
if value>2.5*mean(refv)
    value=mean(refv);
end
         end
         end
         refv=[value refv];
        [pks,locs] = findpeaks(y,'MINPEAKDISTANCE',round(minsep*Fs),'MINPEAKHEIGHT',coef*value);
        locs=locs+j;
%     peaks=[peaks pks];
%     locations=[locations locs];

% %Falsepeakelimination
%    figure() 
%        subplot(2,1,1)
%         plot(j+1:k,y)
%         hold on
%         scatter(locs,pks,'m')
%         xlabel('Sample Number')
%         ylabel({'(a)','d(k)'})
%         axis tight
%         set(gca,'ylim', [0 1.5],'yTick',[0:0.5:1.5],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:5000:k],'fontsize',12)
        
 refdiff=[];
for a=2:length(locs)
difference=locs(a)-locs(a-1);
refdiff=[refdiff difference];
 
end 
counterloc=[refdiff counterloc];
refloc=mean(refdiff);
refstd=[std(counterloc) refstd];
std(refstd);
if std(counterloc)>obtimum
decoeff=0.5;
else
    decoeff=0.8;
end
 
diffmin=decoeff*refloc;    
        
 
for a=1:length(refdiff)
    if a+5<=length(refdiff)
        
         if pks(a+1)>prep*value
             refdiff(a)=refdiff(a);
         else
    if refdiff(a)<diffmin && refdiff(a+1)<diffmin
            alter=refdiff(a)+refdiff(a+1);
            if alter>diffmin
     
                refdiff(a)=0;
                refdiff(a+1)=alter;
                a=a+1;
            else
                alter=refdiff(a)+refdiff(a+1)+refdiff(a+2);
                if alter>diffmin
                   
                refdiff(a)=0;
                refdiff(a+1)=0;
                refdiff(a+2)=alter;
                a=a+2;
                else
                alter=refdiff(a)+refdiff(a+1)+refdiff(a+2)+refdiff(a+3);
                if alter>diffmin
                   
                refdiff(a)=0;
                refdiff(a+1)=0;
                refdiff(a+2)=0;
                refdiff(a+3)=alter;
                a=a+3;
                else
                    alter=refdiff(a)+refdiff(a+1)+refdiff(a+2)+refdiff(a+3)+refdiff(a+4);
                                    if alter>diffmin
                   
                refdiff(a)=0;
                refdiff(a+1)=0;
                refdiff(a+2)=0;
                refdiff(a+3)=0;
                refdiff(a+4)=alter;
                a=a+4;
                                    else
                          alter=refdiff(a)+refdiff(a+1)+refdiff(a+2)+refdiff(a+3)+refdiff(a+4)+refdiff(a+5);
                                    if alter>diffmin
                   
                refdiff(a)=0;
                refdiff(a+1)=0;
                refdiff(a+2)=0;
                refdiff(a+3)=0;
                refdiff(a+4)=0
                refdiff(a+5)=alter;
                a=a+5;
                                    end
                                        
                                    end
                end
                end
            end
 
        end
    end
    end
    end
 
refdiffex=[];
refpeka1=[];
for a=1:length(refdiff)
if refdiff(a)>0
    refdiffex=[refdiffex refdiff(a)];
else
end
 
end
peka=[];
peka(1)=pks(1);
for a=1:length(pks)-1
    if refdiff(a)>0
        peka=[peka pks(a+1)];
    else
    end
end
loca=[];
loca(1)=locs(1);
 
for b=1:length(refdiffex)
 loca(b+1)=loca(b)+refdiffex(b);
end
 
%        subplot(2,1,2)
%         plot(j+1:k,y)
%         hold on
%         scatter(loca,peka,'m')
%         xlabel('Sample Number')
%         ylabel({'(b)','w(k)'})
%         axis tight
%         set(gca,'ylim', [0 1.5],'yTick',[0:0.5:1.5],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:5000:k],'fontsize',12)
%   %endFalsepeakelimination      
        j=k;
 
    finalpeaks=[finalpeaks peka];
    finallocations=[finallocations loca];
end
%% False Peak elimination due to overlapping
 
b=diff(finallocations);
ap=mean(b);
locax=[];
peakx=[];
excludedloc=[];
excludedpeak=[];
stand_dev=std(b);
variance=var(b);
v=mean(refv);
 
I=[];
 
for i=1:length(finallocations)
    if i>1 && i<length(finallocations)
        if finalpeaks(i) >1.5*v
            locax=[locax finallocations(i)];
            peakx=[peakx finalpeaks(i)]; 
        else
            differearl=b(i-1);
           differ=b(i);
        if differ<=aperture*ap
            
            if finalpeaks(i)<finalpeaks(i+1)
            excludedloc=[excludedloc finallocations(i)];
            excludedpeak=[excludedpeak finalpeaks(i)];
            I=[i I];
            elseif differearl<=differ
                 excludedloc=[excludedloc finallocations(i)];
            excludedpeak=[excludedpeak finalpeaks(i)];
            else
            excludedloc=[excludedloc finallocations(i+1)];
            excludedpeak=[excludedpeak finalpeaks(i+1)];
            I=[i I];
            end
        else
        locax=[locax finallocations(i)];
        peakx=[peakx finalpeaks(i)];
        end
        end
    else
        locax=[locax finallocations(i)];
        peakx=[peakx finalpeaks(i)];
    end
        
end
     Beatsbeforesearchback=length(locax)
        
        exed=excludedloc;
        exep=excludedpeak;
        finaldifer=diff(locax)';
        I;
  %%   Search Back   
 pdmean=mean(finaldifer);
 pdstd=std(finaldifer);
 differarray=[];
 
 storepk=[];
 storeloc=[];
 difr2=[];
 difr1=[];
 tictoc=[];
 
    coef2=0.1;
    coef3=0.75;
 
 
 if pdstd<100
for i=2:length(locax)
       dfr1=locax(i)-locax(i-1);
       
      if dfr1>2*round(pdmean)-4*round(pdstd) 
          y=Sb(locax(i-1):locax(i));
          [sbpks,sblocs] = findpeaks(y,'MINPEAKDISTANCE',abs(round(pdmean)-2*round(pdstd)));
              
          for j=1:length(sblocs)   
              sblocs(j)=sblocs(j)+locax(i-1);
              tictoc=[tictoc sblocs(j)];
              difr2(j)=locax(i)-sblocs(j);
              difr1(j)=sblocs(j)-locax(i-1);
              if sbpks(j)>coef2*mean(refv)
                  if difr2(j)>round(pdmean)*coef3 && difr1(j)>round(pdmean)*coef3
          storepk=[storepk sbpks(j)];
          storeloc=[storeloc sblocs(j)];
 
              end
          end
      end
end
end
                
 sblocax=[locax storeloc];
 sbpeakx=[peakx storepk];
 
 if length(sblocax)>length(locax)
 %bb=length(sblocax)-length(locax);
 sqlocax=sort(sblocax);
 sqpeakx=[];
 
for qq=1:length(sqlocax)
    sqpeakx=[sqpeakx Sb(sqlocax(qq))];
        
end
 
 j=0;
sectloc=[];
sectpeak=[];
op=1;
count=0;
n=floor(length(SS)/N);
 for i=1:n
sectloc=[];
sectpeak=[];
 
             k=i*N;
        y=Sb(j+1:k);
        z=SS(j+1:k);
        z1=Sa(j+1:k);
        for m=op:length(sqlocax)
        if sqlocax(m)<k
            sectloc=[sectloc sqlocax(m)];
            sectpeak=[sectpeak sqpeakx(m)];
            
        elseif sqlocax(m)>k
            op=m;
            break
            
        else
        end
        end
        
%     figure()
%         subplot(3,1,1)
%         plot(j+1:k,z1)
          %titile('(a)')
%         axis tight
%         subplot(3,1,2)
%         plot(j+1:k,z)
          %titile('(b)')
%         axis tight
%        subplot(3,1,3)
         %titile('(c)')
%        axis tight
%         plot(j+1:k,y)
%         hold on
%         scatter(sectloc,sectpeak,'m')
%         xlabel('(c) Missed Peaks Detected')
          
%         ylabel('Amplitude')
%         axis tight
 
        j=k;
        count=count+length(sectloc);
        
 end
 Totalbeats=length(sqlocax)
                                else
 Totalbeats=length(locax)
 
 j=0;
sectloc=[];
countloc=[];
sectpeak=[];
countpeak=[];
op=1;
count=0;
n=floor(length(SS)/N);
 for i=1:n
sectloc=[];
sectpeak=[];
 
             k=i*N;
        y=Sb(j+1:k);
        z=SS(j+1:k);
        z1=Sa(j+1:k);
        for m=op:length(locax)
        if locax(m)<k
            sectloc=[sectloc locax(m)];
            sectpeak=[sectpeak peakx(m)];
        elseif locax(m)>k
            op=m;
            break 
        else
        end
        end
 
        
%     figure()
%         subplot(3,1,1)
%         plot(j+1:k,z1)
%         axis tight
%         subplot(3,1,2)
%         plot(j+1:k,z)
%         axis tight
%        subplot(3,1,3)
%         plot(j+1:k,y)
%         hold on
%         scatter(sectloc,sectpeak,'m')
%         xlabel('(c) Missed Peaks Detected')
%         ylabel('Amplitude')
%         axis tight
        j=k;
        count=count+length(sectloc); 
 end
 end
 else
     Totalbeats=length(locax);
 end
 toc;

 
 %% Final placement
 
if pdstd<100
    if length(sblocax)>length(locax)
    solemn=sqlocax;
    solemn2=sqpeakx;
    else
        solemn=locax;
        solemn2=peakx;
    end
else
    solemn=locax;
    solemn2=peakx;
end
 j=0;
op=1;
mp=1;
count=0;
N=25000;
n=floor(length(SS)/N);
 for i=1:n
rgh=[];
rgh2=[];
fl=[];
pl=[];
 
             k=i*N;
             y=Sa(j+1:k);
             y1=SS(j+1:k);
        z=Sc(j+1:k);
        z1=Sb(j+1:k);
        for m=op:length(solemn)
        if solemn(m)<k
            rgh=[rgh solemn(m)];
           rgh2=[rgh2 solemn2(m)];
        elseif solemn(m)>k
            op=m;
            break 
        else
        end
        end
        for m=mp:length(finallocations)
                if finallocations(m)<k
            fl=[fl finallocations(m)];
          pl=[pl finalpeaks(m)];
        elseif finallocations(m)>k
            mp=m;
            break 
        else
        end
        end
 maxvalue=[];
 minvalue=[];
 for pp=1:length(rgh)
     if rgh(pp)-100>0 && rgh(pp)+100<length(SS)
     [maxvalue(pp) maxindex(pp)]=max(SS((rgh(pp)-100):(rgh(pp)+100)));
%     if -1*Sa==gior1
 %    else
     [minvalue(pp) maxindex(pp)]=min(SS((rgh(pp)-100):(rgh(pp)+100)));
     
     if abs(minvalue(pp))>2.5*maxvalue(pp)
         maxvalue(pp)=minvalue(pp);
     end
  %   end
     elseif rgh(pp)+100>length(SS)
          [maxvalue(pp) maxindex(pp)]=max(SS((rgh(pp)-100):(rgh(pp))));
     elseif rgh(pp)-100<0
          [maxvalue(pp) maxindex(pp)]=max(SS((rgh(pp)):(rgh(pp)+100)));
     end
 end

%  %Search back       
%    figure()
%         subplot(2,1,1)
%         plot(j+1:k,z1)
%         hold on
%         scatter(fl,pl,'m')
%         xlabel('Sample Number')
%         ylabel({'(a)','y(k)'})
%         axis tight
%                 set(gca,'ylim', [-0.2 1.6],'yTick',[-0.2:0.8:1.6],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:12500:k],'fontsize',12)
%         subplot(2,1,2)
%         plot(j+1:k,z1)
%         hold on
%         scatter(rgh,rgh2,'m')
%         xlabel('Sample Number')
%         ylabel({'(b)','z(k)'})
%         axis tight
%         set(gca,'ylim', [-0.2 1.6],'yTick',[-0.2:0.8:1.6],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:12500:k],'fontsize',12)
% %Search back

% %preprocessing
%         figure()
%         subplot(3,1,1)
%         plot(j+1:k,y)
%         xlabel('Sample Number')
%         %title('(a)')
%          ylabel({'(a)','x(k)'})
%         axis tight
%         set(gca,'ylim', [-0.6 0.6],'yTick',[-0.6:0.4:0.6],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:500:k],'fontsize',12)
%         
%         subplot(3,1,2)
%         plot(j+1:k,y1)
%         xlabel('Sample Number')
%          ylabel({'(b)','x_2(k)'})
%         axis tight
%         set(gca,'ylim', [-0.5 1],'yTick',[-0.5:0.5:1],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:500:k],'fontsize',12)
%         
%         subplot(3,1,3)
%         plot(j+1:k,z1)
%         xlabel('Sample Number')
%          ylabel({'(c)',['x_3(k)']})
%         axis tight
%        set(gca,'ylim', [-0.1 0.4],'yTick',[-0.1:0.2:0.4],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:500:k],'fontsize',12)
%    %endpreprocessing     

% %peakdetection
%         figure()
%          subplot(3,1,1)
%          plot(j+1:k,y)
%          xlabel('Sample Number')
%          ylabel({'(a)','x(k)'})
%         axis tight
%         set(gca,'ylim', [-0.6 0.6],'yTick',[-0.6:0.4:0.6],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:500:k],'fontsize',12)
%          subplot(3,1,2)
%         plot(j+1:k,z1)
%         hold on
%         scatter(rgh,rgh2,'r')
%         xlabel('Sample Number')
%          ylabel({'(b)','x_3(k)'})
%         axis tight
%         set(gca,'ylim', [-0.1 0.4],'yTick',[-0.1:0.2:0.4],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:500:k],'fontsize',12)
%         subplot(3,1,3)
%         plot(j+1:k,y)
%         hold on
%         stem(rgh,maxvalue,'r')
%         xlabel('Sample Number')
%          ylabel({'(c)','z(k)'})
%        axis tight
%         set(gca,'ylim', [-0.6 1],'yTick',[-0.6:0.6:1],'fontsize',12)
%         set(gca,'xlim', [j k],'xTick',[j:500:k],'fontsize',12)
%  %endpeakdetection
 
%% Final results
        figure()
        subplot(3,1,1)
        plot(j+1:k,z1)
        hold on
        scatter(rgh,rgh2,'m')
        xlabel('Sample Number')
        
         ylabel({'(a)','x_3(k)'})
        axis tight
        %set(gca,'ylim', [-0.2 1],'yTick',[-0.2:0.5:1],'fontsize',12)
        %set(gca,'xlim', [j k],'xTick',[j:2500:k],'fontsize',12)
        subplot(3,1,3)
        plot(j+1:k,y)
        hold on
        stem(rgh,maxvalue)
        xlabel('Sample Number')
         ylabel({'(c)','x(k)'})
         
        axis tight
        %set(gca,'ylim', [-1.5 2.5],'yTick',[-1.5:1.5:2.5],'fontsize',12)
        %set(gca,'xlim', [j k],'xTick',[j:2500:k],'fontsize',12)
       subplot(3,1,2)
        plot(j+1:k,y1)
        hold on
        scatter(rgh,maxvalue,'r')
        xlabel('Sample Number')
       % xlabel('(c) Missed Peaks Detected')
        ylabel({'(b)','x_2(k)'})
        
        axis tight
        %set(gca,'ylim', [-1.5 2.5],'yTick',[-1.5:1.5:2.5],'fontsize',12)
        %set(gca,'xlim', [j k],'xTick',[j:2500:k],'fontsize',12)
%Finalresults
        
        j=k;
        count=count+length(rgh); 
 end
