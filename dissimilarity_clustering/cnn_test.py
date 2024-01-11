# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from keras.models import Sequential
from keras.layers import Dense, BatchNormalization, Dropout
from keras.optimizers import Adam, SGD, RMSprop, Adadelta, Adagrad, Adamax, Nadam, Ftrl
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.wrappers.scikit_learn import KerasClassifier
from math import floor
from sklearn.metrics import make_scorer, accuracy_score
from bayes_opt import BayesianOptimization
from sklearn.model_selection import StratifiedKFold
from keras.layers import LeakyReLU
LeakyReLU = LeakyReLU(alpha=0.1)
import warnings
import os
warnings.filterwarnings('ignore')
pd.set_option("display.max_columns", None)

# Make scorer accuracy
score_acc = make_scorer(accuracy_score)

# Specify the directory containing your files
directory_path = "all_pd"

# Initialize the master list
master_list = []
labels_list = []  # List to store labels

# Iterate through each file in the directory
for i, filename in enumerate(os.listdir(directory_path)):
    if filename.endswith('.hist'):
        # Construct the full file path
        file_path = os.path.join(directory_path, filename)

        # Load the data from the file
        data = np.loadtxt(file_path)

        # Extract the second column into a list and append to the master list
        second_column = data[:60, 1].tolist()
        master_list.append(second_column)

        # Assign labels based on the condition
        label = (i // 25) % 13 + 1
        labels_list.append(label)

# Convert the master list to a NumPy array
master_array = np.array(master_list)

# Create a DataFrame from master_array
df = pd.DataFrame(data=master_array,
                  index=[i for i in range(master_array.shape[0])],
                  columns=['f'+str(i) for i in range(master_array.shape[1])])

df['label'] = labels_list  # Add labels as a new column

# Separate features (X) and target variable (y)
X = df.drop('label', axis=1)
y = df['label']

# Split the data into training and test sets
test_size = 0.2  # Adjust as needed
random_seed = 42  # Set a random seed for reproducibility
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_seed)

# Create function
def nn_cl_bo(neurons, activation, optimizer, learning_rate,  batch_size, epochs ):
    optimizerL = ['SGD', 'Adam', 'RMSprop', 'Adadelta', 'Adagrad', 'Adamax', 'Nadam', 'Ftrl','SGD']
    optimizerD= {'Adam':Adam(lr=learning_rate), 'SGD':SGD(lr=learning_rate),
                 'RMSprop':RMSprop(lr=learning_rate), 'Adadelta':Adadelta(lr=learning_rate),
                 'Adagrad':Adagrad(lr=learning_rate), 'Adamax':Adamax(lr=learning_rate),
                 'Nadam':Nadam(lr=learning_rate), 'Ftrl':Ftrl(lr=learning_rate)}
    activationL = ['relu', 'sigmoid', 'softplus', 'softsign', 'tanh', 'selu',
                   'elu', 'exponential', LeakyReLU,'relu']
    neurons = round(neurons)
    activation = activationL[round(activation)]
    batch_size = round(batch_size)
    epochs = round(epochs)
    def nn_cl_fun():
        opt = Adam(lr = learning_rate)
        nn = Sequential()
        nn.add(Dense(neurons, input_dim=60, activation=activation))
        nn.add(Dense(neurons, activation=activation))
        nn.add(Dense(1, activation='sigmoid'))
        nn.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
        return nn
    es = EarlyStopping(monitor='accuracy', mode='max', verbose=0, patience=20)
    nn = KerasClassifier(build_fn=nn_cl_fun, epochs=epochs, batch_size=batch_size,
                         verbose=0)
    kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)
    score = cross_val_score(nn, X_train, y_train, scoring=score_acc, cv=kfold, fit_params={'callbacks':[es]}).mean()
    return score

# Set paramaters
params_nn ={
    'neurons': (10, 100),
    'activation':(0, 9),
    'optimizer':(0,7),
    'learning_rate':(0.01, 1),
    'batch_size':(200, 1000),
    'epochs':(20, 100)
}

# Run Bayesian Optimization
nn_bo = BayesianOptimization(nn_cl_bo, params_nn, random_state=111)
nn_bo.maximize(init_points=25, n_iter=4)

params_nn_ = nn_bo.max['params']
activationL = ['relu', 'sigmoid', 'softplus', 'softsign', 'tanh', 'selu',
               'elu', 'exponential', LeakyReLU,'relu']
params_nn_['activation'] = activationL[round(params_nn_['activation'])]
params_nn_

# # Create function
# def nn_cl_bo2(neurons, activation, optimizer, learning_rate, batch_size, epochs,
#               layers1, layers2, normalization, dropout, dropout_rate):
#     optimizerL = ['SGD', 'Adam', 'RMSprop', 'Adadelta', 'Adagrad', 'Adamax', 'Nadam', 'Ftrl','SGD']
#     optimizerD= {'Adam':Adam(lr=learning_rate), 'SGD':SGD(lr=learning_rate),
#                  'RMSprop':RMSprop(lr=learning_rate), 'Adadelta':Adadelta(lr=learning_rate),
#                  'Adagrad':Adagrad(lr=learning_rate), 'Adamax':Adamax(lr=learning_rate),
#                  'Nadam':Nadam(lr=learning_rate), 'Ftrl':Ftrl(lr=learning_rate)}
#     activationL = ['relu', 'sigmoid', 'softplus', 'softsign', 'tanh', 'selu',
#                    'elu', 'exponential', LeakyReLU,'relu']
#     neurons = round(neurons)
#     activation = activationL[round(activation)]
#     optimizer = optimizerD[optimizerL[round(optimizer)]]
#     batch_size = round(batch_size)
#     epochs = round(epochs)
#     layers1 = round(layers1)
#     layers2 = round(layers2)
#     def nn_cl_fun():
#         nn = Sequential()
#         nn.add(Dense(neurons, input_dim=60, activation=activation))
#         if normalization > 0.5:
#             nn.add(BatchNormalization())
#         for i in range(layers1):
#             nn.add(Dense(neurons, activation=activation))
#         if dropout > 0.5:
#             nn.add(Dropout(dropout_rate, seed=123))
#         for i in range(layers2):
#             nn.add(Dense(neurons, activation=activation))
#         nn.add(Dense(1, activation='sigmoid'))
#         nn.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['accuracy'])
#         return nn
#     es = EarlyStopping(monitor='accuracy', mode='max', verbose=0, patience=20)
#     nn = KerasClassifier(build_fn=nn_cl_fun, epochs=epochs, batch_size=batch_size, verbose=0)
#     kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)
#     score = cross_val_score(nn, X_train, y_train, scoring=score_acc, cv=kfold, fit_params={'callbacks':[es]}).mean()
#     return score

# params_nn2 ={
#     'neurons': (10, 100),
#     'activation':(0, 9),
#     'optimizer':(0,7),
#     'learning_rate':(0.01, 1),
#     'batch_size':(200, 1000),
#     'epochs':(20, 100),
#     'layers1':(1,3),
#     'layers2':(1,3),
#     'normalization':(0,1),
#     'dropout':(0,1),
#     'dropout_rate':(0,0.3)
# }

# # Run Bayesian Optimization
# nn_bo = BayesianOptimization(nn_cl_bo2, params_nn2, random_state=111)
# nn_bo.maximize(init_points=25, n_iter=4)

# params_nn_ = nn_bo.max['params']
# learning_rate = params_nn_['learning_rate']
# activationL = ['relu', 'sigmoid', 'softplus', 'softsign', 'tanh', 'selu',
#                'elu', 'exponential', LeakyReLU,'relu']
# params_nn_['activation'] = activationL[round(params_nn_['activation'])]
# params_nn_['batch_size'] = round(params_nn_['batch_size'])
# params_nn_['epochs'] = round(params_nn_['epochs'])
# params_nn_['layers1'] = round(params_nn_['layers1'])
# params_nn_['layers2'] = round(params_nn_['layers2'])
# params_nn_['neurons'] = round(params_nn_['neurons'])
# optimizerL = ['Adam', 'SGD', 'RMSprop', 'Adadelta', 'Adagrad', 'Adamax', 'Nadam', 'Ftrl','Adam']
# optimizerD= {'Adam':Adam(lr=learning_rate), 'SGD':SGD(lr=learning_rate),
#              'RMSprop':RMSprop(lr=learning_rate), 'Adadelta':Adadelta(lr=learning_rate),
#              'Adagrad':Adagrad(lr=learning_rate), 'Adamax':Adamax(lr=learning_rate),
#              'Nadam':Nadam(lr=learning_rate), 'Ftrl':Ftrl(lr=learning_rate)}
# params_nn_['optimizer'] = optimizerD[optimizerL[round(params_nn_['optimizer'])]]
# params_nn_

# def nn_cl_fun():
#     nn = Sequential()
#     nn.add(Dense(params_nn_['neurons'], input_dim=10, activation=params_nn_['activation']))
#     if params_nn_['normalization'] > 0.5:
#         nn.add(BatchNormalization())
#     for i in range(params_nn_['layers1']):
#         nn.add(Dense(params_nn_['neurons'], activation=params_nn_['activation']))
#     if params_nn_['dropout'] > 0.5:
#         nn.add(Dropout(params_nn_['dropout_rate'], seed=123))
#     for i in range(params_nn_['layers2']):
#         nn.add(Dense(params_nn_['neurons'], activation=params_nn_['activation']))
#         nn.add(Dense(1, activation='sigmoid'))
#         nn.compile(loss='binary_crossentropy', optimizer=params_nn_['optimizer'], metrics=['accuracy'])
#     return nn
# es = EarlyStopping(monitor='accuracy', mode='max', verbose=0, patience=20)
# nn = KerasClassifier(build_fn=nn_cl_fun, epochs=params_nn_['epochs'], batch_size=params_nn_['batch_size'],
#                         verbose=0)
# nn.fit(X_train, y_train, validation_data=(X_test, y_test), verbose=1)