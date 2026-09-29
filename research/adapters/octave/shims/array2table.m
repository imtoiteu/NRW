function t = array2table(A, varargin)
% Minimal Octave stand-in for MATLAB's array2table (Octave has no `table`
% class). Returns a struct with one column-vector field per variable name.
% Used only so that unmodified scatter-opt code can run under Octave.
  names = {};
  for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'VariableNames')
      names = varargin{k+1};
    end
  end
  if isempty(names)
    for j = 1:columns(A)
      names{j} = sprintf('Var%d', j);
    end
  end
  t = struct();
  for j = 1:numel(names)
    t.(names{j}) = A(:, j);
  end
end
