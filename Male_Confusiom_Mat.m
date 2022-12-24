clc; clear all; close all

dirc_cr_path = '.\Refrences\Child';
dirc_fr_path = '.\Refrences\Female';
dirc_mr_path = '.\Refrences\Male';

dirc_cr = dir(fullfile(dirc_cr_path,'*R.wav'));
dirc_fr = dir(fullfile(dirc_fr_path,'*R.wav'));
dirc_mr = dir(fullfile(dirc_mr_path,'*R.wav'));

dirc_ct_path = '.\Train\Male';
dirc_ct = dir(fullfile(dirc_ct_path,'*.wav'));

%%%%% Determine if the user is C/F/M based on (alhamdullilah) word%%%%%

%%%%% child reference %%%%%%
F1 = fullfile(dirc_cr_path,dirc_cr(end).name);     
[cr,fs1] = audioread(F1);    
if(fs1 ~= 16000)
    audiowrite(F1,cr,16000);
end    

[coeffs_cr,delta,deltaDelta] = mfcc(cr,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_cr_mat = [coeffs_cr, delta, deltaDelta];
coeffs_cr_mat = coeffs_cr_mat - sum(coeffs_cr_mat)/size(coeffs_cr_mat,1);
coeffs_cr_mat(~isfinite(coeffs_cr_mat))=0;
coeffs_cr_mat = fillmissing(coeffs_cr_mat,'linear',2,'EndValues','nearest'); 

%%%%% female reference %%%%%%
F2 = fullfile(dirc_fr_path,dirc_fr(end).name);     
[fr,fs2] = audioread(F2);    
if(fs2 ~= 16000)
    audiowrite(F2,fr,16000);
end 

[coeffs_fr,delta,deltaDelta] = mfcc(fr,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_fr_mat = [coeffs_fr, delta, deltaDelta];
coeffs_fr_mat = coeffs_fr_mat - sum(coeffs_fr_mat)/size(coeffs_fr_mat,1);
coeffs_fr_mat(~isfinite(coeffs_fr_mat))=0;
coeffs_fr_mat = fillmissing(coeffs_fr_mat,'linear',2,'EndValues','nearest'); 

%%%%% male reference %%%%%%
F3 = fullfile(dirc_cr_path,dirc_cr(end).name);     
[mr,fs3] = audioread(F3);    
if(fs3 ~= 16000)
    audiowrite(F3,mr,16000);
end    
[coeffs_mr,delta,deltaDelta] = mfcc(mr,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_mr_mat = [coeffs_mr, delta, deltaDelta];
coeffs_mr_mat = coeffs_mr_mat - sum(coeffs_mr_mat)/size(coeffs_mr_mat,1);
coeffs_mr_mat(~isfinite(coeffs_mr_mat))=0;
coeffs_mr_mat = fillmissing(coeffs_mr_mat,'linear',2,'EndValues','nearest'); 

%%%%% choosing a reference based on (alhamdullilah) word%%%%%
tot_dist = [];
Correct_ref_mapping = 0;
ref_list = [];
user_list = [];
for i = 1:47:numel(dirc_ct)
dirc_ct_user = dirc_ct(i:i+46);
F = fullfile(dirc_ct_path,dirc_ct_user(end).name);   
[ct,fs] = audioread(F);    
if(fs ~= 16000)
    audiowrite(F,ct,16000);        
end 

[coeffs_ct,delta,deltaDelta] = mfcc(ct,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_ct_mat = [coeffs_ct, delta, deltaDelta];     
coeffs_ct_mat = coeffs_ct_mat - sum(coeffs_ct_mat)/size(coeffs_ct_mat,1);
coeffs_ct_mat(~isfinite(coeffs_ct_mat))=0;
coeffs_ct_mat = fillmissing(coeffs_ct_mat,'linear',2,'EndValues','nearest'); 

len = max(size(coeffs_cr_mat,1), size(coeffs_ct_mat,1));
coeffs_ct_mat(end+1:len,:) = 0;
coeffs_cr_mat(end+1:len,:) = 0;

[dtw_dist, ix, iy] = dtw(coeffs_cr_mat, coeffs_ct_mat, 'squared');
dist_cr_ct = min(dist(coeffs_cr_mat(ix),coeffs_ct_mat(iy)));

len = max(size(coeffs_fr_mat,1), size(coeffs_ct_mat,1));
coeffs_ct_mat(end+1:len,:) = 0;
coeffs_fr_mat(end+1:len,:) = 0;

[dtw_dist, ix, iy] = dtw(coeffs_fr_mat, coeffs_ct_mat, 'squared');
dist_fr_ct = min(dist(coeffs_fr_mat(ix),coeffs_ct_mat(iy)));

len = max(size(coeffs_mr_mat,1), size(coeffs_ct_mat,1));
coeffs_ct_mat(end+1:len,:) = 0;
coeffs_mr_mat(end+1:len,:) = 0;

[dtw_dist, ix, iy] = dtw(coeffs_mr_mat, coeffs_ct_mat, 'squared');
dist_mr_ct = min(dist(coeffs_mr_mat(ix),coeffs_ct_mat(iy)));

[min_ref, ind] = min([dist_cr_ct,dist_fr_ct,dist_mr_ct]);

if(ind == 1)
    selected_ref_path = dirc_cr_path;
    selected_ref = dirc_cr;
    %disp('The user is mapped to a Child refrence');    
    ref_list = [ref_list 1];
    Correct_ref_mapping = Correct_ref_mapping +1;
elseif (ind == 2)
    selected_ref_path = dirc_fr_path;
    selected_ref = dirc_fr;   
    %disp('The user is mapped to a Female refrence');
    ref_list = [ref_list 2];
    
elseif (ind == 3)
    selected_ref_path = dirc_mr_path;
    selected_ref = dirc_mr;
    %disp('The user is mapped to a Male refrence');
    ref_list = [ref_list 3];
     Correct_ref_mapping = Correct_ref_mapping +1;
end

dists = [];
for k = 1:2:numel(dirc_ct_user)-1

    [refk1, fs] = audioread(fullfile(selected_ref_path,selected_ref(k).name));
    if(fs ~= 16000)
        audiowrite(fullfile(selected_ref_path,selected_ref(k).name),refk1,16000);        
    end 
    [refk1, fs] = audioread(fullfile(selected_ref_path,selected_ref(k).name));

    
    [refk2, fs] = audioread(fullfile(selected_ref_path,selected_ref(k+1).name));
    if(fs ~= 16000)
        audiowrite(fullfile(selected_ref_path,selected_ref(k+1).name),refk2,16000);        
    end 
    [refk2, fs] = audioread(fullfile(selected_ref_path,selected_ref(k+1).name));

    
    [testk1, fs] = audioread(fullfile(dirc_ct_path,dirc_ct_user(k).name));
    if(fs ~= 16000)
        audiowrite(fullfile(dirc_ct_path,dirc_ct_user(k).name),testk1,16000);        
    end 
    [testk1, fs] = audioread(fullfile(dirc_ct_path,dirc_ct_user(k).name));

    
    [testk2, fs] = audioread(fullfile(dirc_ct_path,dirc_ct_user(k+1).name));
    if(fs ~= 16000)
        audiowrite(fullfile(dirc_ct_path,dirc_ct_user(k+1).name),testk2,16000);        
    end 
    [testk2, fs] = audioread(fullfile(dirc_ct_path,dirc_ct_user(k+1).name));

    
    %%%%%%%% Word 1 from Test %%%%%%%%%%%%
    s1 = refk1;
    [coeffs,delta,deltaDelta] = mfcc(s1,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    coeffs_mat = [coeffs, delta, deltaDelta];     
    s1 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

    s2 = testk1;
    [coeffs,delta,deltaDelta] = mfcc(s2,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    coeffs_mat = [coeffs, delta, deltaDelta];     
    s2 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);
     
     s1(~isfinite(s1))=0;
     s2(~isfinite(s2))=0;
     s1 = fillmissing(s1,'linear',2,'EndValues','nearest'); 
     s2 = fillmissing(s2,'linear',2,'EndValues','nearest'); 
    
    len = max(size(s1,1), size(s2,1));
    s1(end+1:len,:) = 0;
    s2(end+1:len,:) = 0;

     [dtw_dist, ix, iy] = dtw(s1, s2, 'squared');
    dist_refk1_testk1 = min(dist(s1(ix),s2(iy)));
    
    s1 = refk2;
    [coeffs,delta,deltaDelta] = mfcc(s1,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    coeffs_mat = [coeffs, delta, deltaDelta];     
    s1 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

    s2 = testk1;
    [coeffs,delta,deltaDelta] = mfcc(s2,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    coeffs_mat = [coeffs, delta, deltaDelta];     
    s2 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);
     
     s1(~isfinite(s1))=0;
     s2(~isfinite(s2))=0;
     s1 = fillmissing(s1,'linear',2,'EndValues','nearest'); 
     s2 = fillmissing(s2,'linear',2,'EndValues','nearest');
     
     len = max(size(s1,1), size(s2,1));
    s1(end+1:len,:) = 0;
    s2(end+1:len,:) = 0;
    
    [dtw_dist, ix, iy] = dtw(s1, s2, 'squared');
    dist_refk2_testk1 = min(dist(s1(ix),s2(iy)));
    
     %%%%%%%% Word 2 from Test %%%%%%%%%%%%
    s1 = refk1;
    [coeffs,delta,deltaDelta] = mfcc(s1,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    s1 = [coeffs, delta, deltaDelta];     
    s1 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

    s2 = testk2;
    [coeffs,delta,deltaDelta] = mfcc(s2,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    s2 = [coeffs, delta, deltaDelta];     
     
     s1(~isfinite(s1))=0;
     s2(~isfinite(s2))=0;
     s1 = fillmissing(s1,'linear',2,'EndValues','nearest'); 
     s2 = fillmissing(s2,'linear',2,'EndValues','nearest'); 
     %s2(isfinite(s2)|isnan(s2)) = 0;
     
     len = max(size(s1,1), size(s2,1));
    s1(end+1:len,:) = 0;
    s2(end+1:len,:) = 0;
    
     [dtw_dist, ix, iy] = dtw(s1, s2,'squared');
    dist_refk1_testk2 = min(dist(s1(ix),s2(iy)));
    
    s1 = refk2;
    [coeffs,delta,deltaDelta] = mfcc(s1,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    coeffs_mat = [coeffs, delta, deltaDelta];     
    s1 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

    s2 = testk2;
    [coeffs,delta,deltaDelta] = mfcc(s2,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
    coeffs_mat = [coeffs, delta, deltaDelta];     
    s2 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);
    
     s1(~isfinite(s1))=0;
     s2(~isfinite(s2))=0;
     s1 = fillmissing(s1,'linear',2,'EndValues','nearest'); 
     s2 = fillmissing(s2,'linear',2,'EndValues','nearest'); 
     
     len = max(size(s1,1), size(s2,1));
    s1(end+1:len,:) = 0;
    s2(end+1:len,:) = 0;
    
    [dtw_dist, ix, iy] = dtw(s1, s2, 'squared');
    dist_refk2_testk2 = min(dist(s1(ix),s2(iy)));
    
    %%%%%%%%%%% 4 Distances for each pair %%%%%%%%%%%%%%
    dists = [dists dist_refk1_testk1 dist_refk2_testk1 dist_refk2_testk2 dist_refk1_testk2];
    
end
tot_dist = [tot_dist; dists];
end

conf_mat_ct = [];
Missmaches =[];
user_num = [];
%tot_Missmaches =[];
THRESHOLD = sum(tot_dist,1)/size(tot_dist,1);
for i = 1: 4: 92
    correct_word11 = 0;
    wrong_word12 = 0;
    wrong_word21 = 0;
    correct_word22 = 0;
    wrong_other1 = 0;
    wrong_other2 = 0;
    
for j = 1:size(tot_dist,1)
    if tot_dist(j,i) < tot_dist(j,i+1)
        correct_word11 = correct_word11 + 1;
    elseif tot_dist(j,i) > tot_dist(j,i+1) && tot_dist(j,i) < THRESHOLD(i)
        wrong_word12 = wrong_word12 + 1;  
        Missmaches = [Missmaches ceil(i/2)]; % word index in user's data
        user_num = [user_num j];
    else wrong_other1 = wrong_other1 + 1;
    end
    
   if tot_dist(j,i+2) < tot_dist(j,i+3)
       correct_word22 = correct_word22 + 1;        
    elseif tot_dist(j,i+2) > tot_dist(j,i+3) && tot_dist(j,i+2) < THRESHOLD(i+1)
        wrong_word21 = wrong_word21 + 1;
        Missmaches = [Missmaches ceil(i/2)];
        user_num = [user_num j];
    else wrong_other2 = wrong_other2 + 1;
    end
end
%tot_Missmaches = [tot_Missmaches; j Missmaches];
conf_mat_ct = [conf_mat_ct; ceil(i/4) correct_word11 wrong_word12 wrong_other1 wrong_word21 correct_word22 wrong_other2];
end
conf_mat_ct

Pair = conf_mat_ct(:,1);
Word1 = ones([size(conf_mat_ct,1),1]);
correct_Word11 = conf_mat_ct(:,2);
wrong_Word12 = conf_mat_ct(:,3);
wrong_other1 = conf_mat_ct(:,4);
Word2 = ones([size(conf_mat_ct,1),1])*2;
wrong_Word21 = conf_mat_ct(:,5);
correct_Word22 = conf_mat_ct(:,6);
wrong_other2 = conf_mat_ct(:,7);

correct_score = sum(correct_Word11)+sum(correct_Word11);
total_words = sum(wrong_Word12)+sum(wrong_Word21)+sum(correct_Word11)+sum(correct_Word22)...
    +sum(wrong_other2)+sum(wrong_other1);
Acc_Word_mapping = (correct_score/total_words)*100
Acc_ref_mapping = (Correct_ref_mapping/size(tot_dist,1))*100

t = table(Pair,Word1,correct_Word11,wrong_Word12,wrong_other1, ...
  Pair,Word2,wrong_Word21,correct_Word22,wrong_other2);    
uit.Data = t;
uit.ColumnName = {'Pair','Word1','correct_Word11', ...
          'wrong_Word12','wrong_other1','Pair', 'Word2', 'wrong_Word21', ...
          'correct_Word22','wrong_other2'};
uit.RowName = 'numbered';

fig = uifigure;
tbl = uitable(fig);
tbl.Data = uit.Data;

%%%%%%%%%%%%%%% 3 examples of mismatched words %%%%%%%%%%%%%%%
x = [180 50 100];
for i = 1:length(x)
if mod(Missmaches(x(i)),2)==0
ind = Missmaches(x(i));   

ref_num= ref_list(user_num(x(i)));
if ref_num == 1
    selected_ref_path = dirc_cr_path;
    selected_ref = dirc_cr;
elseif ref_num == 2
    selected_ref_path = dirc_fr_path;
    selected_ref = dirc_fr;
elseif ref_num == 3
    selected_ref_path = dirc_mr_path;
    selected_ref = dirc_mr;
end
[ref, fs] = audioread(fullfile(selected_ref_path,selected_ref(ind-1).name));
soundsc(ref,fs);
pause(2.5)  
[test, fs] = audioread(fullfile(dirc_ct_path,dirc_ct((user_num(x(i))-1)*47+ind).name));
soundsc(test,fs);
pause(2.5)  

s1 = ref(n:n+duration-1);
[coeffs,delta,deltaDelta] = mfcc(s1,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_mat = [coeffs, delta, deltaDelta];     
s1 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

s2 = test(n:n+duration-1);
[coeffs,delta,deltaDelta] = mfcc(s2,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_mat = [coeffs, delta, deltaDelta];     
s2 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

 s1(~isfinite(s1))=0;
 s2(~isfinite(s2))=0;
 s1 = fillmissing(s1,'linear',2,'EndValues','nearest'); 
 s2 = fillmissing(s2,'linear',2,'EndValues','nearest'); 
 
len = max(size(s1,1), size(s2,1));
    s1(end+1:len,:) = 0;
    s2(end+1:len,:) = 0;
    
 [dtw_dist, ix, iy] = dtw(s1, s2, 'squared');
 Local_dist = dist(s1(ix),s2(iy));
 t = 1:length(Local_dist);
 figure(i)
 h = plot(t,Local_dist,'LineWidth',1.0);
 hold on
 ld = sprintf('Mismatch (Word 2 is mapped to word 1) in Pair #%d', ind/2);
 xlabel('time')
 ylabel('Local Distance')
 title('Train Data')
 legend(h, ld);         %update the legends
 drawnow(); 
elseif mod(Missmaches(x(i)),2)==1
ind = Missmaches(x(i));   

ref_num= ref_list(user_num(x(i)));
if ref_num == 1
    selected_ref_path = dirc_cr_path;
    selected_ref = dirc_cr;
elseif ref_num == 2
    selected_ref_path = dirc_fr_path;
    selected_ref = dirc_fr;
elseif ref_num == 3
    selected_ref_path = dirc_mr_path;
    selected_ref = dirc_mr;
end
[ref, fs] = audioread(fullfile(selected_ref_path,selected_ref(ind+1).name));
soundsc(ref,fs);
pause(2.5)  
[test, fs] = audioread(fullfile(dirc_ct_path,dirc_ct((user_num(x(i))-1)*47+ind).name));
soundsc(test,fs);
pause(2.5)  

s1 = ref;
[coeffs,delta,deltaDelta] = mfcc(s1,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_mat = [coeffs, delta, deltaDelta];     
s1 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

s2 = test;
[coeffs,delta,deltaDelta] = mfcc(s2,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_mat = [coeffs, delta, deltaDelta];     
s2 = coeffs_mat - sum(coeffs_mat)/size(coeffs_mat,1);

 s1(~isfinite(s1))=0;
 s2(~isfinite(s2))=0;
 s1 = fillmissing(s1,'linear',2,'EndValues','nearest'); 
 s2 = fillmissing(s2,'linear',2,'EndValues','nearest'); 
 
 len = max(size(s1,1), size(s2,1));
    s1(end+1:len,:) = 0;
    s2(end+1:len,:) = 0;
    
 [dtw_dist, ix, iy] = dtw(s1, s2, 'squared');
 Local_dist = dist(s1(ix),s2(iy));
 t = 1:length(Local_dist);
 figure(i)
 h = plot(t,Local_dist,'LineWidth',1.0);
 hold on
 ld = sprintf('Mismatch (Word 1 is mapped to word 2) in Pair #%d', (ind+1)/2);
 xlabel('time')
 ylabel('Local Distance')
 title('Train Data')
 legend(h, ld);         %update the legends
 drawnow(); 
end
end

