function [filename info]=name_file(info)
% Generator of the filename from structure info (input) following
% agreement. It returns filename (char):
% 'experimentname_notes_codewords_YYMMDD.mat'
% It is returning updated structure info with actual date as well, but must
% be followed by updating of the structure info outside of the function

filename=info.experiment; % stores original name of the experiment


if (isfield(info,'section')) && (strcmpi(info.modality,'mri')) % it`s for MRI only
    filename=[filename '_sec' num2str(info.section,'%02u')];
end

try lengthofnotes = length(info.notes{1}); catch end;
if (exist('lengthofnotes','var'))
    if lengthofnotes==0
      info.notes={};
      clear lengthofnotes;
    end;
end

for n=1:length(info.notes)  % gets and stores notes
    filename=[filename '_' info.notes{n}];
end

for n=1:length(info.codewords)  % gets and stores codewords
    filename=[filename '_' info.codewords{n}];
end

string=datestr(now,'yymmdd');
info.date=string;
filename=[filename '_' info.date '.mat']; % constructs whole filename