function varargout = precond_GUI(varargin)
% PRECOND_GUI M-file for precond_GUI.fig
%      PRECOND_GUI, by itself, creates a new PRECOND_GUI or raises the existing
%      singleton*.
%
%      H = PRECOND_GUI returns the handle to a new PRECOND_GUI or the handle to
%      the existing singleton*.
%
%      PRECOND_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRECOND_GUI.M with the given input arguments.
%
%      PRECOND_GUI('Property','Value',...) creates a new PRECOND_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before precond_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to precond_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help precond_GUI

% Last Modified by GUIDE v2.5 12-Sep-2012 15:22:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @precond_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @precond_GUI_OutputFcn, ...
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


% --- Executes just before precond_GUI is made visible.
function precond_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to precond_GUI (see VARARGIN)

% Choose default command line output for precond_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes precond_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global maxiter

maxiter = 10000;

tol = 10^(-12);
set(handles.tolerance,'string',tol);

set(handles.text14,'backgroundcolor', 0.8*[0.5625    0.8    0.4])
set(handles.text15,'backgroundcolor', [ 0    0.4    1.0000])

axis off

% A = get_matrix(handles);

%Test DA CANCELLARE
% A = randi(3,10);
% A = A'*A;
%%%

A = eye(10);
set(handles.matrixA,'string','eye(10)');



b = randi(5,10,1);
set(handles.vectorb,'string','randi(5,10,1)');


P = eye(10);
set(handles.precondP,'string','eye(10)');

M = P'*A*P;

cla
hold on
c = [0.5625    0.8    0.4];
stairplot(A,c,'no',handles);
c =  [ 0    0.4    1.0000];
stairplot(M,c,'yes',handles);
hold off

update_knumbers(handles);

% --- Outputs from this function are returned to the command line.
function varargout = precond_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function matrixA_Callback(hObject, eventdata, handles)
% hObject    handle to matrixA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of matrixA as text
%        str2double(get(hObject,'String')) returns contents of matrixA as a double

update_knumbers(handles);


% --- Executes during object creation, after setting all properties.
function matrixA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matrixA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function precondP_Callback(hObject, eventdata, handles)
% hObject    handle to precondP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of precondP as text
%        str2double(get(hObject,'String')) returns contents of precondP as a double

update_knumbers(handles);


% --- Executes during object creation, after setting all properties.
function precondP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to precondP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function[] = stairplot(A,c,shift,handles)
hold on
axes(handles.axes1)
%c is the colour of the plot (string of the form eg c = 'r)' 
[n, ~] = size(A);
d = sortrows(eig(A));
if strcmp(shift,'no')
    x = 1:n;
else
    x = 1:n;
    x = x + 0.3;
end
y = d';

%Round eigenvalues to 3 d.p. for stylistic purposes
N = 10^(-2);
lab = round(y/N) * N;
h = bar(handles.axes1,x,y,0.3);
set(h,'FaceColor',c)
axis image
axis square
set(get(h,'BaseLine'),'LineWidth',2,'LineStyle',':')

% --- Executes on selection change in preconditioner_menu.
function preconditioner_menu_Callback(hObject, eventdata, handles)
% hObject    handle to preconditioner_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns preconditioner_menu contents checkas cell array
%        contents{get(hObject,'Value')} returns selected item from preconditioner_menu

update_knumbers(handles);



% --- Executes during object creation, after setting all properties.
function preconditioner_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to preconditioner_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function[x_k, counter] = precconjgradmeth(A,b,P,tol)
global maxiter
%Presets
A = P'*A*P;

size(P);
size(b);

b = P'*b;
x_k = 0;
r_k = b;
d_k = r_k;
accept = 1;
counter = 0;

while and(accept == 1, counter < maxiter);
    counter = counter + 1;
    %Step 2
    v_k = A*d_k;
    w_k = ((norm(r_k))^2/(d_k'*v_k));

    %Step 3
    x_k = x_k + w_k*d_k;
    temp = r_k;
    r_k = r_k - w_k*v_k;

    %Step 4
    if norm(r_k) <= tol;
        accept = 0;
        x_k = P*x_k;
        break
    end

    %Step 5
    b_k = ((norm(r_k))^2)/((norm(temp))^2);
    d_k = r_k + b_k*d_k;

    %Step 6
    %Do nothing!
end

if counter == maxiter
    waitfor(msgbox('WARNING: Maximum number of iterations reached.',...
        'Warning','warning','modal'));
end
    



function vectorb_Callback(hObject, eventdata, handles)
% hObject    handle to vectorb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vectorb as text
%        str2double(get(hObject,'String')) returns contents of vectorb as a double

% --- Executes during object creation, after setting all properties.
function vectorb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vectorb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showbutton.
function showbutton_Callback(hObject, eventdata, handles)
% hObject    handle to showbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    b = get_vectorb(handles);
    tol = get_tolerance(handles);
    [A, P] = get_matrices(handles);
    assert(~isnan(tol) && isreal(tol) && ~isinf(tol) && tol > eps,...
        'Invalid tolerance! (Is it too small?)');
    assert( ndims(b)==2 && ndims(A) == 2 && ndims(P) == 2 &&...
        size(A,1) == size(A,2) && size(A,2) == size(P,1) &&...
        size(P,1) == size(P,2) && size(P,2) == size(b,1) && ...
        size(b,2)==1 && size(P,2)>=2,...
    ['A and P must be square matrices of equal size, '...
    'and b must be a column vector' ...
    'conformable with A and P for multiplication!']);
        
    n = size(A);

    [sol1, counter1] = precconjgradmeth(A,b,eye(n),tol);
    [sol2, counter2] = precconjgradmeth(A,b,P,tol);

    set(handles.noiterA,'string',counter1)
    set(handles.noiterPAP,'string',counter2)

    fprintf('\n --------------------------------------- \n ')
    fprintf(' The exact solution x of your system is: \n')
    x = A\b;
    disp(x);
    fprintf(' The l2 error of the not preconditioned system is: \n')
    disp(norm(x - sol1));
    fprintf(' The l2 error of the preconditioned system is: \n')
    disp(norm(x - sol2));

    M = P'*A*P;


    cla
    hold on
    c = [0.5625    0.8    0.4];
    stairplot(A,c,'no',handles);
    c =  [ 0    0.4    1.0000];
    stairplot(M,c,'yes',handles);
    hold off

    d1 = sortrows(eig(A));
    set(handles.listboxA,'string',num2str(d1));
    d2 = sortrows(eig(M));
    set(handles.listboxPAP,'string',num2str(d2));
catch err
    waitfor(msgbox(err.message,'Error','error','modal'));
end





% --- Executes when selected object is changed in preconditioner_panel.
function preconditioner_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in preconditioner_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

update_knumbers(handles);


% --- Executes on selection change in listboxA.
function listboxA_Callback(hObject, eventdata, handles)
% hObject    handle to listboxA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxA contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxA


% --- Executes during object creation, after setting all properties.
function listboxA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxPAP.
function listboxPAP_Callback(hObject, eventdata, handles)
% hObject    handle to listboxPAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxPAP contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxPAP


% --- Executes during object creation, after setting all properties.
function listboxPAP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxPAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in closebutton.
function closebutton_Callback(hObject, eventdata, handles)
% hObject    handle to closebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);

function tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tolerance as text
%        str2double(get(hObject,'String')) returns contents of tolerance as
%        a double

% --- Executes during object creation, after setting all properties.
function tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_knumbers (handles)
try
    [A, P] = get_matrices(handles);
    assert( ndims(A) == 2 && ndims(P) == 2 &&...
        size(A,1) == size(A,2) && size(A,2) == size(P,1) &&...
        size(P,1) == size(P,2) && size(P,2)>=2,...
    'A and P must be square matrices of equal size!');
    set(handles.condPAP, 'string', num2str(precondition_number(P'*A*P)));
    set(handles.condA, 'string', num2str(precondition_number(A)));
    set(handles.showbutton,'Enable','on');
catch err
    set(handles.showbutton,'Enable','off');
    waitfor(msgbox(err.message,'Error','error','modal'));
end
    

function b = get_vectorb (handles)
b = str2num(get(handles.vectorb, 'string'));

function k = precondition_number ( M )
e = abs(eig(M));
k = max(e)/min(e);


function pre_method = get_preconditioner_method ( handles )
if get(handles.custom_radiobutton, 'value') == 1
    pre_method = 'custom';
else
    pre_method = 'preset';
end

function [A, P] = get_matrices (handles)
A = str2num(get(handles.matrixA,'string'));
A = A'*A;
switch get_preconditioner_method(handles)
    case 'custom'
        P = str2num(get(handles.precondP, 'String'));
    case 'preset'
        [n, ~] = size(A);
        h = get(handles.preconditioner_menu, 'value');
        switch h
            case 1
                P = eye(n);
            case 2 %LU factorisation
                [L,~,~] = lu(A);
                P = inv(L)';
            case 3 %Cholesky factorisation
                S = diag(diag(A));
                Q = chol(S);
                P = inv(Q);
        end
end

function tol = get_tolerance ( handles )
tol = str2double(get(handles.tolerance, 'String'));
