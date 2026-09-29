function h = height(t)
% Minimal Octave stand-in for MATLAB's height() (not implemented in Octave 8),
% applied to the struct that replaces a MATLAB table in the scatter-opt adapter:
% number of rows of the first field.
  if isstruct(t)
    fn = fieldnames(t);
    h = rows(t.(fn{1}));
  else
    h = rows(t);
  end
end
