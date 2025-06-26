function genJpattern(fnc, x)
     
% --------------------------------------------------------------------------------------   
% generate and save new JPattern for arbitrary input and state
% fnc inputs
    % fnc - function handle, ODE right hand side rhs(x)
    % x   - state vector


% -------------------------------------------------------------------------------------- 

% Jacobian calculation
nstates = length(x);
jac = sparse([],[],[],nstates,nstates);
for i1=1:nstates
    dx = sparse(i1,1,1e-99*1i,nstates,1);
    jac(:,i1) = (imag(fnc( x+dx ))*1e99)~=0;
end

% save JPattern to MATLAB function myJPattern()
[ijac,jjac]=find(jac);
fid = fopen(which('myJPattern'),'w');
fprintf(fid,'function out = myJPattern()\n');
fprintf(fid,'A=['); 
fprintf(fid,'%u %u; ',[ijac,jjac]');
fprintf(fid,'];\n');
fprintf(fid,'out = sparse(A(:,1),A(:,2),(ones([%u,1])),%u,%u);',numel(ijac),nstates,nstates);
fclose(fid);

end
