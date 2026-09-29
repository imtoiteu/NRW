function import(varargin)
% No-op stand-in for MATLAB's `import` (not implemented in Octave 8; calling it
% raises an error). scatter-opt calls `import scatter_opt.*` inside every
% function; the adapter instead puts the +scatter_opt folder itself on the
% path so the unqualified calls resolve. Nothing else is changed.
end
