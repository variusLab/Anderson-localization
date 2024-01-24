## Localisation d'Anderson d'atomes froids dans un potentiel lumineux désordonné

On considère un système d'atomes froids sans intéraction en présence d'un potentiel lumineux désordonné (champ de tavelures laser invariant dans le temps) et d'un couplage spin-orbite de type Rashba. 

On s'intéresse à la fonction d'onde $\psi(r,t)$ d'une particule qui contient toute l'information sur la dynamique du système.
À deux dimnsions, son évolution est régie par l'équation de Schrödinger (adimensionnée) :
$$i \partial_t \psi(\mathbf{r},t) = \left[ \mathbf{k}^2 + V(\mathbf{r}) + \lambda_R \left(k_y \sigma_x - k_x \sigma_y \right) \right]\psi(\mathbf{r},t)$$
où $\sigma_x$, $\sigma_y$ sont les matrices de Pauli, $\lambda_R > 0$ est l'intensité d'interaction spin-orbite et $\psi = (\psi_+,\psi_-)$ est un spineur.
Les détailles de construction du potentiel désordonné $V(\mathbf{r})$ sont présentés dans le rapport.

Le programme résout cette équation en utilisant la méthode _split-step Fourier symétrique_. 
Il produit 4 figures:

- figure 1 (animation): évolution dans le temps de la densité de probabilité de présence du système ($|\psi(x,y,t)|^2$);
- figure 2: capture de la figure 1 (profile 2D de $|\psi(x,y,t)|^2$) à un instant $t=t_{max}/2$;
- figure 3: évolution temporelle de la taille moyenne du paquet d'onde $\ell(t) = \sqrt{( \langle\psi(r,t)|\mathbf{r}^2|\psi(r,t)\rangle}$);
- figure 4: évolution temporelle de la probabilité de présence du système dans tout l'espace (on vérifie que $\langle\psi|\psi\rangle(t)=1$ pour tout $t$).

<img src = "M 128-dt 0.005-tmax_7.0-sigma 5.0-V0 0.1-lambdaR_3.0.gif">


