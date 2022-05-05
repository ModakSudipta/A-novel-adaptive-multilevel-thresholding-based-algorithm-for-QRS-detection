close all;
clear all;clc;
x=rdsamp('104.dat');
x1=x(:,1);
x2=x(:,2);
flag=1;
fs=360;
Sa=x1;
%% preprocessing
Sa=Sa./max(Sa);

    SS = highpass(Sa,5,fs);
    SS= lowpass(SS,35,fs);
      %[SS baseline]=removebaseline(approx,fs,0.1);
      
       for i=1:length(SS)
     if SS(i)<0
         SS(i)=0.8*SS(i);
     else
         SS(i)=SS(i);
     end
       end
 Sb=abs(SS);
 Sb=smooth(Sb,round(0.05*fs));
 %Sb=Sb./max(Sb);
 %Sb=Sb.^2;
 
 
 %% Segmentation
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
else
end
 
coef1=0.4;
aperture=0.48;
%% Peak Detection

      [pks,locs] = findpeaks(Sb,'MINPEAKDISTANCE',round(0.32*fs));
      
 peak=[];
loca =[];
SIG_LEV = 0; 
diff_RR =[]; 
mean_RR = 0;
sum=0;
avger=mean(pks);
THR_SIG =avger  ;
THR_NOISE =0.6*avger;
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;
tresloc=[];
d=diff(locs);
beater=[];
COUNT=0;
COUNT2=0;
[pks,locs] = findpeaks(Sb,'MINPEAKDISTANCE',round(0.3*fs),'MINPEAKHEIGHT',0.6*avger);
%avger2=mean(pks);
cheater=[];
for i=1:length(locs)


    if pks(i)>=THR_SIG 
        if length(loca)>5
            amp=pks(i);
                peak=[peak amp];
                place=locs(i);
                loca= [loca place];
                SIG_LEV=0.1*amp+0.9*SIG_LEV;
                beatrate=(fs*60)/(loca(end)-loca(end-1));
                beater=[beater beatrate];
                
        else
            amp=pks(i);
                peak=[peak amp];
                place=locs(i);
                loca= [loca place];
                SIG_LEV=0.1*amp+0.9*SIG_LEV;
        end
 elseif pks(i)>=THR_NOISE && pks(i)<THR_SIG
                if length(loca)>5
                     br=(fs*60)/(locs(i)-loca(end));
                 if br<1.5*mean(beater(3:end)) 
                        
                amp=pks(i);
                peak=[peak amp];
                place=locs(i);
                loca= [loca place];
                SIG_LEV=0.1*amp+0.9*SIG_LEV;
                beatrate=(fs*60)/(loca(end)-loca(end-1));
                beater=[beater beatrate];
                 else
                      noisepk(i)=pks(i);
                    noiseloc(i)= locs(i);
                    NOISE_LEV=0.2*noisepk(i)+0.8*NOISE_LEV;   
                 end
                else
                 amp=pks(i);
                peak=[peak amp];
                place=locs(i);
                loca= [loca place];
                SIG_LEV=0.1*amp+0.9*SIG_LEV;
                end
           
         end
         
        
         THR_SIG = 0.6*NOISE_LEV + 0.4*(SIG_LEV - NOISE_LEV);
         THR_NOISE=0.5*THR_SIG;
        
         tres1y(i)=THR_SIG;
         tres1x(i)=locs(i);
         tres2y(i)=THR_NOISE;
         tres2x(i)=locs(i);
         tresloc=[tresloc locs(i)];
          
        
        
end
        Total_beatsafterpd=length(loca)
          
%% False Peak elimination due to high frequency noise
n=floor(length(SS)/N);
j=0;
op=1;
locations=[];
peaks=[];
finalpeaks=[];
finallocations=[];
refv=[];
counterloc=[];
refstd=[];
for i=1:n
 rgh=[];
 rgh2=[];
        k=i*N;
        y=Sb(j+1:k);
        z=Sa(j+1:k);
 for m=op:length(loca)
        if loca(m)<k
            rgh=[rgh loca(m)];
           rgh2=[rgh2 peak(m)];
        elseif loca(m)>k         
            op=m;
            break 
        else
        end
 end
        
 refdiff=diff(rgh);
 counterloc=[refdiff counterloc];
refloc=mean(refdiff);
refstd=[std(refdiff) refstd];
dior=std(refstd);
if std(counterloc)>obtimum
decoeff=0.5;
else
    decoeff=0.7;
end
 %decoeff=0.5;
diffmin=decoeff*refloc;    
 threshavg=mean(rgh2);       
 
for a=1:length(refdiff)
    
    if a+2<=length(refdiff)
        if rgh2(a+1)>0.8*threshavg
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
pekb=[];
pekb(1)=rgh2(1);
for a=1:length(rgh2)-1
    if refdiff(a)>0
        pekb=[pekb rgh2(a+1)];
    else
    end
end
locb=[];
locb(1)=rgh(1);
 
for b=1:length(refdiffex)
 locb(b+1)=locb(b)+refdiffex(b);
end  
      
    finalpeaks=[finalpeaks pekb];
    finallocations=[finallocations locb];
    
%          figure()
%          subplot(2,1,1)
%          plot(j+1:k,y)
%          hold on
%         scatter(rgh,rgh2,'m')
%         xlabel('Sample Number')
%         ylabel({'(b)','w(k)'})
%         axis tight
%         subplot(2,1,2)
%         plot(j+1:k,y)
%         hold on
%         scatter(locb,pekb,'r')
%         xlabel('Sample Number')
%         ylabel({'(b)','w(k)'})
%         axis tight
        %set(gca,'ylim', [0 1.5],'yTick',[0:0.5:1.5],'fontsize',12)
        %set(gca,'xlim', [j k],'xTick',[j:5000:k],'fontsize',12)
  
          j=k;
end

Total_beatsafter_fpe=length(finallocations)
        
        
%% Searchback        
sbrealpeak=[];
sbrealloc=[];

            diff_RR=diff(finallocations);
            dmean_RR=mean(diff_RR);
            dstand=std(diff_RR);
           if dstand<220 %|| dstand<120 
 j=0;
op=1;

N=25000;
n=floor(length(SS)/N);
 for i=1:n

sbfl=[];
sbpl=[];
 
             k=i*N;
        for m=op:length(finallocations)
                if finallocations(m)<k
            sbfl=[sbfl finallocations(m)];
            sbpl=[sbpl finalpeaks(m)];
            
        elseif finallocations(m)>k
            op=m;
            break 
        else
        end
        end
        diff_rr=diff(sbfl);
        dmean_rr=mean(diff_rr);
        for c=1:length(diff_rr)
             if diff_rr(c)>1.75*round(dmean_rr) %>=2*dmean_RR
                 dsb=Sb((sbfl(c)+0.25*dmean_rr):(sbfl(c+1)-0.5*dmean_rr));
          [sbpks,sblocs] = findpeaks(dsb,'MINPEAKDISTANCE',0.65*dmean_rr);
          
          for k=1:length(sbpks)
              
              if sbpks(k)>0.25*avger
              sblocs(k)=sblocs(k)+sbfl(c)+0.25*dmean_rr;
              sbrealpeak=[sbrealpeak sbpks(k)];
              sbrealloc=[sbrealloc sblocs(k)];
              end
              end
          end
             end
            end
          
           end
        
        Total_beats=length(finallocations)+length(sbrealloc)
        
         figure()
         plot(Sb)
        axis tight
        hold on
        scatter(loca,peak)
        hold on
        plot(tres1x,tres1y,'r')
        hold on
        plot(tres2x,tres2y,'k')
        hold on
        scatter(sbrealloc,sbrealpeak,'k')
        
     %% Finalplacement
     
 j=0;
op=1;
mp=1;
count=0;
N=25000;
n=floor(length(SS)/N);
 for i=1:n

fl=[];
pl=[];
f2=[];
p2=[];
p2l=[];
 
             k=i*N;
             
             y=SS(j+1:k);
             z=Sb(j+1:k);

        for m=mp:length(finallocations)
                if finallocations(m)<k
            fl=[fl finallocations(m)];
            pl=[pl finalpeaks(m)];
            p2l=[p2l 1];
        elseif finallocations(m)>k
            mp=m;
            break 
        else
        end
        end
      for m=op:length(sbrealloc)
                if sbrealloc(m)<k && sbrealloc(m)>=j
            f2=[f2 sbrealloc(m)];
            p2=[p2 sbrealpeak(m)];
        elseif sbrealloc(m)>k
            op=m;
            break 
        else
        end
      end  
      figure()
      subplot(3,1,1)
      plot(j+1:k,y)
       subplot(3,1,2)
        plot(j+1:k,y)
        hold on
        stem(fl,p2l,'b')
        hold on
        stem(f2,p2,'k')
        xlabel('Sample Number')
        ylabel({'(a)','y(k)'})
         set(gca,'ylim', [-0.75 1.25])

        subplot(3,1,3)
        plot(j+1:k,z)
        hold on
        scatter(fl,pl,'m')
        hold on
        scatter(f2,p2,'k')
        xlabel('Sample Number')
        ylabel({'(b)','z(k)'})
     
 j=k;
 end