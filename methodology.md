
From our model, we get start with $e_{+}^*,e_{\times}^*, \sigma^*$, then, for $\omega \sim 0.1$ assuming the error in PSF is of the order of $10\%$, we can define  

$$\Delta e_1 = \omega e_1^{PSF} = e_1^{*}-e_1^{PSF},$$ 

$$\Delta e_2 = \omega e_2^{PSF} = e_2^{*}-e_2^{PSF},$$

which are e.g. radial models (qualitatively observed in previous surveys), and $\sigma^{PSF} \sim \sigma^*$ (as an example, a reasonable $\sigma$ would be 0.4 arcsec (or a FWHM of 0.7 arcsec $-$ [median for LSST](https://www.lsst.org/sites/default/files/docs/137.01_Wolff_LSST_System_8x10.pdf)).

We can get $e_1^*, e_2^*$ using $ \begin{bmatrix} e_1 \\ 
e_2
\end{bmatrix} =  \begin{bmatrix}
cos 2\theta & -sin 2\theta \\
sin 2\theta & cos 2\theta 
\end{bmatrix} \begin{bmatrix} e_{+} \\ 
e_{\times} 
\end{bmatrix}$



then $ e_i^{PSF} = e_i^{*}/(1+\omega) $ where $i=1,2$.

By definition: $ e_1^{k} = \frac{M_{xx}^k - M_{yy}^k}{TrM^k} $ and $ e_2^k = \frac{2M_{xy}^k }{TrM^k} $
where $k$ is either PSF or *; and $TrM = M_{xx} + M_{yy}$.

Solving for $M_{ij}$ we get:

$$M^{k} = \frac{TrM^k}{2} \begin{bmatrix}
     e_1 + 1       &  e_2  \\
     e_2       &  - e_1 + 1
\end{bmatrix}$$

then after taking the arithmetic mean of each matrix element $M_{ij}^k$ to get $\langle M^k \rangle = \sum_{l=1}^N M_{ij,l}^k$ for an exposure $l$ over all exposures at a certain point, we can go back to e-space:

$$ e_1^{k} = \frac{M_{xx}^k - M_{yy}^k}{TrM^k} $$ and $$ e_2^k = \frac{2M_{xy}^k }{TrM^k} $$ and $$\sigma = \sqrt{TrM^k/2} $$
