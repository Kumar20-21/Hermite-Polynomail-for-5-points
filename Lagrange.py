def lagrange(t, n, x, f):  # t is the point of evaluation, n is the number of nodes, x is an array of node ponts and f is the array of function values corresonding to node points.
    f_x = 0.0
    for ii in range(n):
        prod = 1.0
        for jj in range(n):  # This loop will construct l_i(x)
            if jj == ii:  # Skipping the ith factor in the product
                continue
            else:
                prod *= (t-x[jj])/(x[ii]-x[jj])  # (t-x_j)/(x_i-x_j) term in the l_i(x)
        f_x += f[ii]*prod  # Multiplying the value of function at node with l_i(x)
    return f_x
