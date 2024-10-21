# Code adapted from David & Christine
# Summer 2023
import scipy.optimize as opt
import numpy as np
from astropy.cosmology import FlatLambdaCDM


def fit_sigma_mu2(params, z, x_0, errors, cov, pec=250, dint=0.1):
    dmuz = pec / 3e5
    dlens = 0.055 * z
    dmb = -2.5 * 0.434 * errors[:, 1] / x_0  # mb = 2.5log10(x0) + constant
    dx1 = params[0] * errors[:, 2]
    dc = params[1] * errors[:, 3]
    dcx1 = 2 * params[0] * params[1] * cov[:, 2, 3]
    dmbx1 = 2 * params[0] * cov[:, 1, 2] * (-2.5 / (x_0 * np.log(10.0)))
    dmbc = 2 * params[1] * cov[:, 1, 3] * (-2.5 / (x_0 * np.log(10.0)))
    err = np.sqrt(dint ** 2 + dmuz ** 2 + dlens ** 2 + dmb ** 2 + dx1 ** 2 + dc ** 2 + dcx1 + dmbx1 + dmbc)
    return err
def fit_mu(params, z, x_0, x_1, c, masses):
    alpha = params[0]
    beta = params[1]
    gamma = params[2]
    bigM = params[3]
    # bigM = params[3:5]
    m_shift = np.zeros(x_0.shape[0])
    if masses is not None:
        m_shift[masses > 10] = 1

    mu_base = -2.5 * np.log10(x_0) + 10.635 + (alpha * x_1) - (beta * c) + gamma * m_shift
    mu_base -= bigM
    # mu_base[(z > 0.033) & (z <= 0.067)] -= bigM[1]
    # mu_base[(z > 0.067) & (z <= 0.1)] -= bigM[2]
    return mu_base
def cost_mu(params, z, x_0, x_1, c, errors, cov, masses):
    expected = np.array(FlatLambdaCDM(70, 0.3).distmod(z).value)
    fit = np.array(fit_mu(params[:4], z, x_0, x_1, c, masses))
    err = fit_sigma_mu2(params[:4], z, x_0, errors, cov, dint=params[4])
    n = z.shape[0]
    return np.sum(-(fit - expected) ** 2 / (2 * (err ** 2)) - np.log(np.sqrt(2 * np.pi * err ** 2))) * -1
def optimize_alpha_beta(path: str):
    """
    Using SALT3 paramater file, optimizes for alpha & beta.
    :param path: str; location of SALT3 paramater file.
    :return: results: dict; results of optimization function
    """
    # Load data
    data = np.genfromtxt(path, delimiter=', ', dtype=str, skip_header=1)
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]

    # Set arrays
    z = data[:, hdr.index('z')].astype(float)
    x_0 = data[:, hdr.index('x0')].astype(float)
    x_1 = data[:, hdr.index('x1')].astype(float)
    c = data[:, hdr.index('c')].astype(float)
    masses = data[:, hdr.index('hostMass')].astype(float)

    # Errors of x0, x1, c in a nested numpy array
    e = np.transpose(np.array([np.full(len(x_0), 0.00),
                     data[:, hdr.index('x0_err')].astype(float),
                     data[:, hdr.index('x1_err')].astype(float),
                     data[:, hdr.index('c_err')].astype(float)]))

    # # Errors of empty covariances
    # cv = np.full((len(x_0), 4, 4), 0.00)

    # Load covariances
    cv = np.full((len(x_0), 4, 4), 0.00)
    for k in range(len(z)):
        for i in range(0, 4):
            for j in range(0, 4):
                cv[k, i, j] = data[k, hdr.index('cov'+str(i)+str(j))].astype(float)


    # Put through optimization function
    r = opt.minimize(cost_mu, x0=[0.14, 3.0, 19.36, 0.0, 0.1], args=(z, x_0, x_1, c, e, cv, masses))
    results = r.x
    errs = np.sqrt(np.diag(r.hess_inv))

    return {'alpha': results[0], 'beta': results[1]}
