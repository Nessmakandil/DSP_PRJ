function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 20-Dec-2022 14:03:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirc_cr_path = '.\Refrences\Child';
dirc_fr_path = '.\Refrences\Female';
dirc_mr_path = '.\Refrences\Male';

dirc_cr = dir(fullfile(dirc_cr_path,'*R.wav'));
dirc_fr = dir(fullfile(dirc_fr_path,'*R.wav'));
dirc_mr = dir(fullfile(dirc_mr_path,'*R.wav'));

%%%%% child reference %%%%%%
F1 = fullfile(dirc_cr_path,dirc_cr(end).name);     
[cr,fs1] = audioread(F1);    
if(fs1 ~= 16000)
    audiowrite(F1,cr,16000);
end    

[coeffs_cr,delta,deltaDelta] = mfcc(cr,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_cr_mat = [coeffs_cr, delta, deltaDelta];
coeffs_cr_mat = coeffs_cr_mat - sum(coeffs_cr_mat)/size(coeffs_cr_mat,1);

%%%%% female reference %%%%%%
F2 = fullfile(dirc_fr_path,dirc_fr(end).name);     
[fr,fs2] = audioread(F2);    
if(fs2 ~= 16000)
    audiowrite(F2,fr,16000);
end    
[coeffs_fr,delta,deltaDelta] = mfcc(fr,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_fr_mat = [coeffs_fr, delta, deltaDelta];
coeffs_fr_mat = coeffs_fr_mat - sum(coeffs_fr_mat)/size(coeffs_fr_mat,1);
%%%%% male reference %%%%%%
F3 = fullfile(dirc_cr_path,dirc_cr(end).name);     
[mr,fs3] = audioread(F3);    
if(fs3 ~= 16000)
    audiowrite(F3,mr,16000);
end    
[coeffs_mr,delta,deltaDelta] = mfcc(mr,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
coeffs_mr_mat = [coeffs_mr, delta, deltaDelta];
coeffs_mr_mat = coeffs_mr_mat - sum(coeffs_mr_mat)/size(coeffs_mr_mat,1);

%%%%% choosing a reference based on (alhamdullilah) word%%%%%
tot_dist = [];
Correct_ref_mapping = 0;
ref_list = [];
user_list = [];

[filename, pathname] = uigetfile('*.wav', 'Select a wave file','MultiSelect','on');
% Check if the user selected one file or multiple
if iscell(filename)
     nfiles = length(filename);     
     tot_dist = [];
     for i = 1:47:nfiles
        cr_tt_dis =[];
        fr_tt_dis =[];
        mr_tt_dis =[];

        nwavfile = fullfile(pathname, filesep, filename{1,end});
        [tt,fs] = audioread(nwavfile);    
        if(fs ~= 16000)
            audiowrite(nwavfile,tt,16000);        
        end   
        [coeffs_tt,delta,deltaDelta] = mfcc(tt,16000,'NumCoeffs',12,'OverlapLength',0.01*16000,'WindowLength',0.025*16000);
        coeffs_tt_mat = [coeffs_tt, delta, deltaDelta];     
        coeffs_tt_mat = coeffs_tt_mat - sum(coeffs_tt_mat)/size(coeffs_tt_mat,1);
        
        coeffs_cr_mat(~isfinite(coeffs_cr_mat))=0;
        coeffs_tt_mat(~isfinite(coeffs_tt_mat))=0;
        coeffs_cr_mat = fillmissing(coeffs_cr_mat,'linear',2,'EndValues','nearest'); 
        coeffs_tt_mat = fillmissing(coeffs_tt_mat,'linear',2,'EndValues','nearest'); 
        
        len = max(size(coeffs_cr_mat,1), size(coeffs_tt_mat,1));
        coeffs_tt_mat(end+1:len,:) = 0;
        coeffs_cr_mat(end+1:len,:) = 0;

        [dtw_dist, ix, iy] = dtw(coeffs_cr_mat, coeffs_tt_mat, 'squared');
        dist_cr_tt = min(dist(coeffs_cr_mat(ix),coeffs_tt_mat(iy)));

        coeffs_fr_mat(~isfinite(coeffs_fr_mat))=0;
        coeffs_fr_mat = fillmissing(coeffs_fr_mat,'linear',2,'EndValues','nearest'); 
        
        len = max(size(coeffs_fr_mat,1), size(coeffs_tt_mat,1));
        coeffs_tt_mat(end+1:len,:) = 0;
        coeffs_fr_mat(end+1:len,:) = 0;

        [dtw_dist, ix, iy] = dtw(coeffs_fr_mat, coeffs_tt_mat, 'squared');
        dist_fr_tt = min(dist(coeffs_fr_mat(ix),coeffs_tt_mat(iy)));

        coeffs_mr_mat(~isfinite(coeffs_mr_mat))=0;
        coeffs_mr_mat = fillmissing(coeffs_mr_mat,'linear',2,'EndValues','nearest'); 
        
        len = max(size(coeffs_mr_mat,1), size(coeffs_tt_mat,1));
        coeffs_tt_mat(end+1:len,:) = 0;
        coeffs_mr_mat(end+1:len,:) = 0;
        
        [dtw_dist, ix, iy] = dtw(coeffs_mr_mat, coeffs_tt_mat, 'squared');
        dist_mr_tt = min(dist(coeffs_mr_mat(ix),coeffs_tt_mat(iy)));

        [min_ref, ind] = min([dist_cr_tt,dist_fr_tt,dist_mr_tt]);

        
        if(ind == 1)
            selected_ref_path = dirc_cr_path;
            selected_ref = dirc_cr;
            d = sprintf('selected reference is Child');
            set(handles.text,'String',d);
        elseif (ind == 2)
            selected_ref_path = dirc_fr_path;
            selected_ref = dirc_fr;
            d = sprintf('selected reference is Female');
            set(handles.text,'String',d);
        else
            selected_ref_path = dirc_mr_path;
            selected_ref = dirc_mr;
            d = sprintf('selected reference is Child');
            set(handles.text,'String',d);
        end

        dists = [];
        for k = 1:2:nfiles-1

            [refk1, fs] = audioread(fullfile(selected_ref_path,selected_ref(k).name));
            if(fs ~= 16000)
                audiowrite(fullfile(selected_ref_path,selected_ref(k).name),refk1,16000);        
            end 

            [refk2, fs] = audioread(fullfile(selected_ref_path,selected_ref(k+1).name));
            if(fs ~= 16000)
                audiowrite(fullfile(selected_ref_path,selected_ref(k+1).name),refk2,16000);        
            end 

            [testk1, fs] = audioread(fullfile(pathname, filesep, filename{1,k}));
            if(fs ~= 16000)
                audiowrite(fullfile(pathname, filesep, filename{1,k}),testk1,16000);        
            end 

            [testk2, fs] = audioread(fullfile(pathname, filesep, filename{1,k+1}));
            if(fs ~= 16000)
                audiowrite(fullfile(pathname, filesep, filename{1,k+1}),testk2,16000);        
            end 
            
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
            
            s1 = refk1;
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
            dist_refk1_testk2 = min(dist(s1(ix),s2(iy)));
                        
            dists = [dists dist_refk1_testk1 dist_refk2_testk1 dist_refk2_testk2 dist_refk1_testk2];

        end
        tot_dist = [tot_dist; dists]; 
        tot_dist
     end
    conf_mat_tt = [];
    Missmaches =[];
    
    for i =1:4:92
        correct_word11 = 0;
        wrong_word12 = 0;
        wrong_word21 = 0;
        correct_word22 = 0;
        wrong_other1 = 0;
        wrong_other2 = 0;
        THRESHOLD = sum(tot_dist(i:i+3))/4
    for j = 1:size(tot_dist,1)    
    if tot_dist(j,i) < tot_dist(j,i+1)
        correct_word11 = correct_word11 + 1;
    elseif tot_dist(j,i) > tot_dist(j,i+1) && tot_dist(j,i) < THRESHOLD
        wrong_word12 = wrong_word12 + 1;  
        Missmaches = [Missmaches ceil(i/2)]; % word index in user's data
    else wrong_other1 = wrong_other1 + 1;
    end
    
   if tot_dist(j,i+2) < tot_dist(j,i+3)
       correct_word22 = correct_word22 + 1;        
    elseif tot_dist(j,i+2) > tot_dist(j,i+3) && tot_dist(j,i+2) < THRESHOLD
        wrong_word21 = wrong_word21 + 1;
        Missmaches = [Missmaches ceil(i/2)];
    else wrong_other2 = wrong_other2 + 1;
   end
    end
    conf_mat_tt = [conf_mat_tt; ceil(i/4) correct_word11 wrong_word12 wrong_other1 wrong_word21 correct_word22 wrong_other2];
    end
  Missmaches     
set(handles.axes1,'Units','pixels');
axes(handles.axes1);
ax = handles.axes1;
for v = 1:3
    if mod(Missmaches(v),2)==0
    ind = Missmaches(v); 
        [ref, fs] = audioread(fullfile(selected_ref_path,selected_ref(ind-1).name));
    soundsc(ref,fs);
    pause(2.5) 
    [test, fs] = audioread(fullfile(pathname, filesep, filename{1,ind}));
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
     h(v) = plot(ax,t,Local_dist,'LineWidth',1.0);
     hold on
     ld{v} = sprintf('Mismatch (Word 2 is mapped to word 1) in Pair #%d', ind/2);
     xlabel('time')
     ylabel('Local Distance')
     title('Test Data')
     legend(h(1:v), ld{1:v});         %update the legends
     drawnow(); 
    elseif mod(Missmaches(v),2)==1
    ind = Missmaches(v);  
    
    [ref, fs] = audioread(fullfile(selected_ref_path,selected_ref(ind+1).name));
    soundsc(ref,fs);
    pause(2.5)  
    [test, fs] = audioread(fullfile(pathname, filesep, filename{1,ind}));
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
     h(v) = plot(ax,t,Local_dist,'LineWidth',1.0);
     hold on
     ld{v} = sprintf('Mismatch (Word 1 is mapped to word 2) in Pair #%d', (ind+1)/2);
     legend(h(1:v), ld{1:v});         %update the legends
     xlabel('time')
     ylabel('Local Distance')
     title('Test Data')
     drawnow(); 
    end
  end 
    
    Pair = conf_mat_tt(:,1);
    Word1 = ones([size(conf_mat_tt,1),1]);
    correct_Word11 = conf_mat_tt(:,2);
    wrong_Word12 = conf_mat_tt(:,3);
    wrong_other1 = conf_mat_tt(:,4);
    Word2 = ones([size(conf_mat_tt,1),1])*2;
    wrong_Word21 = conf_mat_tt(:,5);
    correct_Word22 = conf_mat_tt(:,6);
    wrong_other2 = conf_mat_tt(:,7);

    correct_score = sum(correct_Word11)+sum(correct_Word11);
    total_words = sum(wrong_Word12)+sum(wrong_Word21)+sum(correct_Word11)+sum(correct_Word22)...
        +sum(wrong_other2)+sum(wrong_other1);
    Acc_Word_mapping = (correct_score/total_words)*100;
    
    d = sprintf('Accuracy of Word mapping %f %',Acc_Word_mapping);
    set(handles.text1,'String',d);
   
    t = table(Pair,Word1,correct_Word11,wrong_Word12,wrong_other1, ...
      Pair,Word2,wrong_Word21,correct_Word22,wrong_other2);    
    uit.Data = t;
    uit.ColumnName = {'Pair','Word1','correct_Word11', ...
              'wrong_Word12','wrong_other1','Pair', 'Word2', 'wrong_Word21', ...
              'correct_Word22','wrong_other2'};
    uit.RowName = 'numbered';
    fig = uifigure;
    hObject = uitable(fig,'Data', uit.Data);
    guidata(hObject, handles);    
    
 
end

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
