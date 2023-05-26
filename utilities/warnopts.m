function unassigned = warnopts(unassigned)
% warnopts - warn about unassigned options.  rem=warnopts(assignopts(...))

if (length(unassigned))
%  warning(['unrecognized options:', sprintf(' %s', unassigned{1:2:end})])
  error(['unrecognized options:', sprintf(' %s', unassigned{1:2:end})])
end
