import matplotlib.pyplot as plt
import numpy as np
import sncosmo

def SALT2Fit():
    # Register sample data
    data = sncosmo.load_example_data()
    print(data)

    # Create Model
    model = sncosmo.Model(source='salt2')

    # Fit data to model
    result, fitted_model = sncosmo.fit_lc(data,
                                          model,
                                          ['z', 't0', 'x0', 'x1', 'c'],  # parameters of model to vary
                                          bounds={'z': (0.3, 0.7)})  # bounds on parameters (if any)

    # Print Results
    print("Number of chi^2 function calls:", result.ncall)
    print("Number of degrees of freedom in fit:", result.ndof)
    print("chi^2 value at minimum:", result.chisq)
    print("model parameters:", result.param_names)
    print("best-fit values:", result.parameters)
    print("The result contains the following attributes:\n", result.keys())

    # Plot data with fit
    sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)

    #
    model.set(z=0.5)  # set the model's redshift.
    result, fitted_model = sncosmo.fit_lc(data, model,
                                          ['t0', 'x0', 'x1', 'c'])
    sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)

    return