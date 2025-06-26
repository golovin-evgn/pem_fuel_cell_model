function myvec = struct2vec(mystruct,variableNames,variableDimensions)
%struct2vec constructs the state vector from the state struct
%  using state information from cell array 'states'

%%

% initialize counter for state variables
mycounter = 0;

% get size of state description
nStates = length(variableNames);

nTime = size(mystruct.(variableNames{1}),2);

% initialize state vector      
myvec = coder.nullcopy( zeros([sum(variableDimensions),nTime],'like',mystruct.(variableNames{1})) );

for i1 = 1:nStates % for all state variables

  % calculate number of elements to be copied to state vector
  n = variableDimensions(i1);

  % copy data from state struct to state vector
  myvec(mycounter+(1:n),1:nTime) = mystruct.(variableNames{i1});

  % increase state variable counter
  mycounter = mycounter + n;

end

end
