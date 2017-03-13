function sec = testpatchcontact(tpfi,tpei,varargin)
% testpatchcontact(tpi) tests contactpatchcheck() and patchcontact() on each
% type in tpi
%
%   INPUTS
%
%   The parameter tpi is the corresponding type of patches that are extracted 
%   from testcontactpatches(). The testpatchcontact() uses tpi=1:15.
%
%   OUTPUTS
%
%   Return sec is the total execution time.
%
%   OPTIONS
%
%   A list of (name,value) pairs may be given to set options.  Unrecognized 
%   names cause warnings.  If a name is given more than once the last-given 
%   (rightmost) value takes precedence. Extra args are passed on to
%   contactpatchcheck() and patchcontact(). The testpatchcontact(tpi) uses
%   default contactpatchcheck() and patchcontact() options.
%
%
%   EXAMPLES
%
%   >> testpatchcontact(1:2,1:15);
%
% Copyright (C) 2016- Dimitrios Kanoulas

% patchplot options
ppo = {'bw',2,'gw',1,'df',0,'da',0,'dg',0}; dpf = 0; dpe = 1; dbg = 1;

% patchcontact options (various options)
pco={};
%pco = {'dorot',0,'drot',0,'dcp',0,'dcpts',0,'dbg',0}; %no debugging for timing
%pco = {'dorot',0,'drot',0,'dcp',1,'dcpts',1,'dbg',1};
%pco = {'dorot',1,'drot',0,'dcp',1,'dcpts',1,'dbg',1};
%pco = {'dorot',1,'drot',0,'dcp',1,'dcpts',1,'dbg',2};

% end of configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2 types of foot patches and 15 types of environment patches
if (nargin < 1); tpfi = 1:2; tpei = 1:15; end

% setup figure window
persistent fig;
if (isempty(findall(fig,'Type','Figure'))); fig = figure(); end
if (ishandle(fig)); clf(fig,'reset'); % make current w/o focus
else; figure(fig);
end
set(fig,'Color','w'); isempty(fig), ishandle(fig) 

  % Test function for the foot and environment patches
  function ok = test(pf, pe)
  % Given a foot patch pf and an environment patch pe: first check for contact
  % validity and then find the exact contact patch.
  if (contactpatchcheck(pf,pe,varargin{:}))
    [pc,phi,cp,ct,sec] = patchcontact(pf,pe,pco{:},varargin{:});
    if (dbg); fprintf('Time for all contacts=%g\n',sec); end
  end
  ok = 1;
  end

% run the test
testcontactpatches(@test,tpfi,tpei,dpf,dpe,ppo);

end