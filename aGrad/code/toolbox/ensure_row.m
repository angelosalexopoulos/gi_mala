function r = ensure_row(v)
    if iscolumn(v), r = v'; else, r = v; end
end