% gui1.m --------------------------------------------------------
%
% Usually plt is used to build gui applications which include plotting,
% however this example shows that you can use it to create a gui without
% any plots. If you are new to Matlab GUIs, this example is a good one
% to start with since it is simple enough to get the basic ideas. The
% gui objects used here are entirely uicontrols, although the idea of
% a pseudo object is also introduced by using the PseudoSlider object.
% (A PseudoSlider is an object defined by plt consisting of a collection
% of 5 uicontrols designed to work together to control a single parameter.)
%
% This GUI doesn't actually perform any useful function other than to
% demonstrate how to create and the various controls and move them around
% until the GUI appears as desired. The slider callback generates new
% random numbers for the listbox, textbox, and uitable. The remaining
% callbacks are just stubs that notify you that you clicked on the object.
%
% You can most easily absorb the point of this example (and the following
% one called gui2.m) by watching my video which you can find here:
%   www.mennen.org\plt\video\MatlabGUIbuilding.avi
%
% Please address any questions or comments you may have about this example,
% the video, or plt in general to me at the email address shown below.
%
% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function gui1()
  figure('name','gui1','menu','none','pos',[60 60 430 350],'color',[0 .1 .2]);
  cho = {'choice A' 'choice B' 'choice C'}; % choices for popup control
  p = {[.020 .920 .300     ];  % PseudoSlider 1
       [.350 .920 .300     ];  % PseudoSlider 4
       [.680 .920 .300     ];  % PseudoSlider 3
       [.020 .500 .480 .280];  % uitable
       [.540 .500 .440 .280];  % frame
       [.680 .710 .170 .050];  % popup
       [.570 .610 .380 .060];  % slider
       [.570 .520 .170 .060];  % button
       [.780 .520 .170 .060];  % checkbox
       [.020 .050 .480 .400];  % listbox (80 lines)
       [.540 .050 .440 .400]}; % text (10 lines)
% First we create the pseudo sliders and the uitable. Even though we don't use these
% handles we save them anyway because in a real gui we would almost always need them.
  h1 = plt('slider',p{1}, 10,'PseudoSlider 1',@CBsli);
  h2 = plt('slider',p{2}, 60,'PseudoSlider 2',@CBsli);
  h3 = plt('slider',p{3},800,'PseudoSlider 3',@CBsli); 
  h4 = uitable('units','norm','pos',p{4});
% Next we create all 7 uicontrols in one line. Note how the properties for each uicontrol
% in the next three lines (style, string, callback) line up directly under the uicontrol
% command that created the object. 
  h =            [uicontrol uicontrol uicontrol uicontrol uicontrol uicontrol uicontrol];
  set(h,{'style'},{'frame' ;'popup';  'slider';'pushb'  ;'checkbox';'listbox';'text';},...
       {'string'},{'frame1'; cho   ;  'slider';'button1';'check001';''       ;''     },...
       { 'callb'},{''      ; @CBpop;  @CBsli  ;@CBpush  ; @CBcheck ;''       ;''     },...
       'backgr',[.5 1 1],'units','norm',{'pos'},p(5:end));
  set(h(1),'backgr',[1 1 2]/6,'foregr',[1 1 1]/2); % colors for the frame
  set(gcf,'user',[h1 h2 h3 h4 h]);  CBsli;     % save the handles for the slider callback.
%end function gui1

function CBpop(a,b)
                      disp('popup callback');
function CBcheck(a,b)
                      disp('checkbox callback');
function CBpush(a,b)
                      disp('pushbutton callback');

function CBsli(a,b)                      % The slider callback -------------------
  h = get(gcf,'user');                   % Get the handle list
  t = 1e20.^(rand(3,80))/1e6;            % Random numbers (with wide dynamic range)
  set(h(10:11),'fontname','courier',...  % Use the same random table of numbers
    'string',prin('3{%6V   }~, ',t));    %  for both the listbox and the textbox
  set(h(4),'data',100*rand(3,2));        % More random numbers for the uitable
