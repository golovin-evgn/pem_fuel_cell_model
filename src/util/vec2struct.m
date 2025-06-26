function varargout = vec2struct(myvec,variableNames_in,variableDimensions_in)
%vec2struct constructs a struct from a composite vector with fields
% according to the individual elements' names and dimensions.
%  In case of vector elements only, the input is allowed to be
%  a matrix with dim(vector) x dim(time).

%%

% check input types
% assert( isnumeric(variableDimensions) , "Input value 'variableDimensions' is not numeric." );
% assert( ismatrix(variableDimensions)  , "Input value 'variableDimensions' is not a matrix." );

% backwards compatibility / single struct output:
% allow array for dimensions and cell of character arrays for names 
if ~iscell(variableDimensions_in)
    variableDimensions = {variableDimensions_in};
else
    variableDimensions = variableDimensions_in;
end
if ischar(variableNames_in{1})
    variableNames = {variableNames_in};
else
    variableNames = variableNames_in;
end

% initialize counter for state variables
mycounter = 0;

% prepare dummy field names
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';

for i0 = 1:nargout

    % get size of state description
    nStates = length(variableNames{i0});

    % construct struct depending on dimension of elements according to variableDimensions
    % e.g. if matrices are allowed, everything is treated as a matrix

    if isvector(variableDimensions{i0}) % for vector

        %get number of time points, array with dim(stateVector) x dim(time)
        nTime = size(myvec,2);

        for i1 = 1:nStates % for all state variables

            % calculate number of elements to be extracted from state vector
            n = variableDimensions{i0}(i1);

            % write variable name and value into state struct
            mystructs.(alphabet(i0)).(variableNames{i0}{i1}) = myvec(mycounter+(1:n),1:nTime);

            % increase state variable counter
            mycounter = mycounter + n;

        end

    elseif numel(variableDimensions{i0})==2*nStates % for matrix

        for i1 = 1:nStates % for all state variables

            % calculate number of elements to be extracted from state vector
            n1 = variableDimensions{i0}(i1,1);
            n2 = variableDimensions{i0}(i1,2);
            ntotal = n1 * n2;

            % write variable name and value into state struct
            mystructs.(alphabet(i0)).(variableNames{i0}{i1}) = reshape(myvec(mycounter+(1:ntotal)),[n1,n2]);

            % increase state variable counter
            mycounter = mycounter + ntotal;

        end

    else % for tensor

        for i1 = 1:nStates % for all state variables

            % calculate number of elements to be extracted from state vector
            n = variableDimensions{i0}(i1,:);
            ntotal = prod(n,2);

            % write variable name and value into state struct
            mystructs.(alphabet(i0)).(variableNames{i0}{i1}) = reshape(myvec(mycounter+(1:ntotal)),n);

            % increase state variable counter
            mycounter = mycounter + ntotal;

        end

    end

end


for i0=1:nargout
    varargout{i0} = mystructs.(alphabet(i0));
end

end

% variableNames={{'a'; 'b'} {'c'}}; variableDimensions={[1;2] [3]}; myvec=[1 2 3 4 5 6].'; [a b] =vec2struct(myvec,variableNames,variableDimensions)
% variableNames={'a'; 'b'}; variableDimensions=[1;2]; myvec=[1 2 3].'; a=vec2struct(myvec,variableNames,variableDimensions)
% variableNames={'a'; 'b'}; variableDimensions=[1;2]; myvec=[1 2 3; 4 5 6].'; a=vec2struct(myvec,variableNames,variableDimensions)