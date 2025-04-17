function v=jUpperTriMatToVec(m,varargin)
% converts the upper-triangular part of a matrix to a vector
%
% IN:
%   m: matrix
%   offset: offset above leading diagonal, fed to triu function
% OUT:
%   v: vector of upper-triangular values
% 
% v1.0 Oct 2009 Jonas Richiardi
% - initial release

switch nargin
    case 1
        offset=1;
    case 2
        offset=varargin{1};
end

% get indices of upper triangular part (Peter Acklam's trick)
[m_i m_j] = find(triu(ones(size(m)), offset));
% copy to vector
v=zeros(numel(m_i),1);
for v_idx=1:numel(m_i)
    v(v_idx)=m(m_i(v_idx),m_j(v_idx));
end