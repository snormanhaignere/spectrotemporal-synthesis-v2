function index = subplot_sub2ind(n_rows, n_cols, row, col)

% returns the appropriate index for a specific row and column of a 2D subplot
index = col + (row-1) * n_cols;