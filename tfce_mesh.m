function tfce = tfce_mesh(faces, t, dh, E, H)
% FORMAT tfce = tfce_mesh(faces, t, dh, E, H)
% Estimate TFCE
% faces     - faces of surface 
% t         - T map 
% dh        - step size (e.g. dh = max(abs(t))/100)
% E         - TFCE parameter for extent
% H         - TFCE parameter for height
%
% Christian Gaser
% $Id: tfce_mesh.m 125 2017-08-23 14:59:44Z gaser $

tfce = zeros(size(t));
t_max = max(abs(t(:)));

n_steps = ceil(t_max/dh);

%-Compute the (reduced) adjacency matrix
%--------------------------------------------------------------------------
A       = spm_mesh_adjacency(faces);
A       = A + speye(size(A));
A1 = A; A2 = A;

clear A

t1 = t; t2 = t;

for j = 1:n_steps

  thresh = j*dh;

  % positive t-values
  T   = (t1 >= thresh);
  ind = find(T);
  T1 = (t >= thresh);
  ind1 = find(T1);
  
  if ~isempty(ind)
    [C, A1] = find_connected_component(A1, T);
    C = C(ind);
    
    for i = 1:max(C)
      M = find(C == i);
      k = length(M);
      tfce(ind1(M)) = tfce(ind1(M)) + power(k,E)*power(thresh,H)*T1(ind1(M));
    end
    
    t1 = t1(ind);  
  end
  
  % negative t-values
  if min(t2(:)) < 0
    T   = (-t2 >= thresh);
    ind = find(T);
    T2   = (-t >= thresh);
    ind2 = find(T2);
  
    if ~isempty(ind)
      [C, A2] = find_connected_component(A2, T);
      C = C(ind);
    
      for i = 1:max(C)
        M = find(C == i);
        k = length(M);
        tfce(ind2(M)) = tfce(ind2(M)) - power(k,E)*power(thresh,H)*T2(ind2(M));
      end
      
    end
    t2 = t2(ind);  
  end

end

%--------------------------------------------------------------------------
function [C, A] = find_connected_component(A, T)
% find connected components 
% FORMAT C = find_connected_component(A,T)
% A        - a [nxn[ (reduced) adjacency matrix
% T        - a [nx1] data vector (using NaNs or logicals), n = #vertices
%
% C        - a [nx1] vector of cluster indices
%
% modified version from spm_mesh_clusters.m 5065 2012-11-16 20:00:21Z guillaume
%


%-Input parameters
%--------------------------------------------------------------------------
if ~islogical(T)
  T   = ~isnan(T);
end
  
A(~T,:) = [];
A(:,~T) = [];

%-And perform Dulmage-Mendelsohn decomposition to find connected components
%--------------------------------------------------------------------------
[p,q,r] = dmperm(A);
N       = diff(r);
CC      = zeros(size(A,1),1);
for i = 1:length(r)-1
  CC(p(r(i):r(i+1)-1)) = i;
end
C       = NaN(numel(T),1);
C(T)    = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
[N,ni]  = sort(N(:), 1, 'descend');
[ni,ni] = sort(ni);
C(T)    = ni(C(T));
