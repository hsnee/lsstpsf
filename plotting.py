def sns_time_series(x_tuple,y_tuple,outputname,errors=0,two=False, *args,**kwargs):
    """seaborn time series, with error-bands"""
    if (type(outputname)==str)|(type(x_tuple)==tuple)|(type(y_tuple)==tuple):
        pass
    else:
        raise TypeError()
    import matplotlib
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns; sns.set_style('darkgrid')
    import seaborn.timeseries
    x, x_label = x_tuple
    y, y_label = y_tuple
    if two==True:
        x2,x_label2 = x_tuple2
        y2,y_label2 = y_tuple2
    def _plot_std_bars(std=None, central_data=None, ci=None, data=None,*args, **kwargs):
        std = errors
        ci = np.asarray((central_data - std, central_data + std))
        kwargs.update({"central_data": central_data, "ci": ci, "data": data})
        seaborn.timeseries._plot_ci_band(*args, **kwargs)
    seaborn.timeseries._plot_std_bars = _plot_std_bars

    plt.figure()
    sns.tsplot(xip,r,err_style='std_bars')
    sns.tsplot(xim,r,err_style='std_bars',color='r')
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylabel(r'$\xi$')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend([r'$\xi_+$',r'$\xi_-$'],bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
    plt.savefig(outputname+'.pdf')
    plt.close()


def whisker_plot(X,Y,U,angles, V=None,C=None, key=(1,False),pixel_scale=1,color_scale=1,units='rad'):
    import angles, matplotlib
    matplotlib.use("pdf")
    import numpy as np, matplotlib.pyplot as plt
    import seaborn as sns;sns.set_style('darkgrid')
    plt.figure()
    key_value, key_name = key
    if V==None:
        V = U.copy()
    if units=='rad':
        pixel_scale,X,Y = np.array(map(lambda x: angles.r2arcs(x),(pixel_scale,X,Y)))
    elif units=='deg':
        pixel_scale,X,Y = np.array(map(lambda x: angles.d2arcs(x),(pixel_scale,X,Y)))
    elif units=='arcs':
        pass
    else:
        raise ValueError('unknown units')
    X, Y = np.array(map(lambda x: x * pixel_scale, (X,Y)))
    std2fwhm = 2.*np.sqrt(2.*np.log(2.))
    Q = plt.quiver(X-np.median(X),Y-np.median(Y),U,V,std2fwhm*C*color_scale,angles=angles, \
                        headlength=0,headaxislength=0,scale=20,cmap='viridis')
    qk = plt.quiverkey(Q, 0.9, 0.9, key_value, r'$|$' + key_name + '$|=$' + str(key_value), \
                        labelpos='E', coordinates='figure')
    plt.title("whisker plot")
    plt.colorbar(label='FWHM')
    #plt.axis('equal')
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    plt.savefig(name+"_whisker.pdf")
    plt.close()

def HSC_style_plots(X,Y,e1a,e2a,e1b,e2b,sigmaa,sigmab):
    import angles, matplotlib
    matplotlib.use("pdf")
    import numpy as np, matplotlib.pyplot as plt, seaborn as sns; sns.set_style('darkgrid')

    X,Y,e1a,e2a,e1b,e2b,sigmaa,sigmab = np.array(map(lambda x:x[1::5],(X,Y,e1a,e2a,e1b,e2b,sigmaa,sigmab)))
    print 'meshing started'
    X,Y = np.meshgrid(X,Y)
    print 'meshing finished'
    X,Y = X-np.mean(X), Y-np.mean(Y)
    X,Y = np.array(map(lambda x: angles.r2arcs(1)*x,(X,Y)))
    X,Y,e1a,e2a,e1b,e2b,sigmaa,sigmab = np.array(X), np.array(Y), np.array(e1a), np.array(e2a), np.array(e1b), np.array(e2b), np.array(sigmaa), np.array(sigmab)
    print np.size(X), np.size(Y), np.size(e1a)

    plt.figure(1)
    plt.scatter(X,Y,c=e1a-e1b,cmap='viridis')
    plt.colorbar()
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    plt.title(r'$\Delta e_1$')
    plt.savefig('Delta_e1.pdf')
    plt.close()

    plt.figure(2)
    plt.scatter(X,Y,c=e2a-e2b,cmap='viridis')
    plt.colorbar()
    plt.title(r'$\Delta e_2$')
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    plt.savefig('Delta_e2.pdf')
    plt.close()

    plt.figure(3)
    plt.scatter(X,Y,c=(sigmaa-sigmab)/(0.5*(sigmaa+sigmab)),cmap='viridis')
    plt.colorbar()
    plt.title(r'$\frac{\sigma_a-\sigma_b}{<\sigma>}$')
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    plt.savefig('Delta_sigma.pdf')
    plt.close()
