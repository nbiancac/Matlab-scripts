%
%  fix_size.m  ver 1.1  July 12, 2012
%
%  by Tom Irvine
%
function[a]=fix_size(a)
sz=size(a);
if(sz(2)>sz(1))
    a=a';
end