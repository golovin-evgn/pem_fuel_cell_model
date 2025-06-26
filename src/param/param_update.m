function [out] = param_update(p,varargin)
%PARAM_UPDATE updates parameter structure
% 
%   param_update(p,p_update) updates fields in p with all fields
%   with matching names in p_update
% 
%   param_update(p,updateFieldNames,updateFieldValues) updates
%   fields in struct p that match names in updateFieldNames with new values
%   from cell array of values fieldValues 

%   non-matching fields in p_update/updateFieldNames are ignored without
%   error message or warning

switch nargin
    case 2
        p_update = varargin{1};

        % copy new values to output structure

        if ~coder.target('MATLAB')

            % copy input structure and update matching fields
            p_fieldnames=fieldnames(p);
            for i1 = 1:numel(p_fieldnames)
                if isfield(p_update,p_fieldnames{i1})
                    p2.(p_fieldnames{i1}) = p_update.(p_fieldnames{i1});
                else
                    p2.(p_fieldnames{i1}) = p.(p_fieldnames{i1});
                end
            end

            out = p2;

        else

            % overwrite new values directly in input structure
            updateFieldNames=fieldnames(p_update);
            for i1 = 1:numel(updateFieldNames)
                if isfield(p,updateFieldNames{i1})
                    p.(updateFieldNames{i1}) = p_update.(updateFieldNames{i1});
                end
            end

            out = p;

        end

    case 3
        updateFieldNames = varargin{1};
        updateFieldValues = varargin{2};

        % copy new values to output structure

        if ~coder.target('MATLAB')

            % copy input structure and update matching fields
            p_fieldnames=fieldnames(p);
            for i1 = 1:numel(p_fieldnames)
                if any(strcmp(updateFieldNames,p_fieldnames{i1}))
                    p2.(p_fieldnames{i1}) = updateFieldValues{i1};
                else
                    p2.(p_fieldnames{i1}) = p.(p_fieldnames{i1});
                end
            end

            out = p2;

        else

            % overwrite new values directly in input structure
            for i1 = 1:numel(updateFieldNames)
                if isfield(p,updateFieldNames{i1})
                    p.(updateFieldNames{i1}) = updateFieldValues{i1};
                end
            end

            out = p;

        end
end


% test.a=1;test.b='a'; test.c=@(x) x.^2; test_update.a=2;test_update.c=@(x) x.^3; test_new=param_update(test,test_update)
% test.a=1;test.b='a'; updateFieldValues={2}; updateFieldNames={'a'}; test_new=param_update(test,updateFieldNames,updateFieldValues)