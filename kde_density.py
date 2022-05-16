def kde_density(adata, x, y, nbins=300, cmap=None, figsize=[4,4]):
    """
    adata: The annotated data matrix.
    x: a gene's expression at x coordinate
    y: a gene's expression at y coordinate
    nbins: number of bins for the kernel density estimation
    
    # usage
    kde_density(adata, x='Gata3',y='Id2',cmap=plt.cm.Greys)
    """
    # libraries
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import kde
    x= np.array(adata[:,x].X.todense()).reshape(-1)
    y= np.array(adata[:,y].X.todense()).reshape(-1)
    
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    plt.figure(figsize=figsize)
    # Make the plot
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto', cmap=cmap)
    plt.show()

    
