import numpy as np
import pickle
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from keras.optimizers import Adam, SGD, RMSprop, Adadelta, Adagrad, Adamax, Nadam, Ftrl
import matplotlib.pyplot as plt


os.chdir("..")

with open("2b_all_pd_equilibrium", 'rb') as pickle_file:
    all_pd = pickle.load(pickle_file)

with open("labels_pd_equilibrium", 'rb') as pickle_file:
    labels = pickle.load(pickle_file)

print(np.shape(all_pd))

stacked_data = np.vstack(all_pd)
print(np.shape(stacked_data))

# Create a DataFrame
all_df = pd.DataFrame(stacked_data)

# Create a DataFrame
all_df['label'] = labels

# Split data into features and labels
X = all_df.drop('label', axis=1)
y = all_df['label']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Standardize the features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Build the neural network model
opt = Adam()
nn = Sequential()
nn.add(Dense(64, input_dim=X_train_scaled.shape[1], activation='softsign'))
nn.add(Dense(32, activation='relu'))
nn.add(Dropout(0.5)) 
nn.add(Dense(32, activation='relu'))
nn.add(Dense(len(np.unique(labels)), activation='softmax'))
nn.compile(loss='sparse_categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

# Train the model
history = nn.fit(X_train_scaled, y_train, epochs=250, batch_size=32, validation_split=0.33)

# Evaluate the model on the test set
loss, accuracy = nn.evaluate(X_test_scaled, y_test)
print(f"Test Accuracy: {accuracy * 100:.2f}%")

plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Mean Squared Error')
plt.legend()
plt.show()


# # Import packages
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from sklearn.model_selection import train_test_split
# from sklearn.model_selection import cross_val_score
# from keras.models import Sequential
# from keras.layers import Dense, BatchNormalization, Dropout
# from keras.optimizers import Adam, SGD, RMSprop, Adadelta, Adagrad, Adamax, Nadam, Ftrl
# from keras.callbacks import EarlyStopping, ModelCheckpoint
# from keras.wrappers.scikit_learn import KerasClassifier
# from math import floor
# from sklearn.metrics import make_scorer, accuracy_score
# from bayes_opt import BayesianOptimization
# from sklearn.model_selection import StratifiedKFold
# from keras.layers import LeakyReLU
# LeakyReLU = LeakyReLU(alpha=0.1)
# import warnings
# warnings.filterwarnings('ignore')
# pd.set_option("display.max_columns", None)

# # Make scorer accuracy
# score_acc = make_scorer(accuracy_score)

# # Create function
# def nn_cl_bo(neurons, activation, optimizer, learning_rate,  batch_size, epochs ):
#     optimizerL = ['SGD', 'Adam', 'RMSprop', 'Adadelta', 'Adagrad', 'Adamax', 'Nadam', 'Ftrl','SGD']
#     optimizerD= {'Adam':Adam(lr=learning_rate), 'SGD':SGD(lr=learning_rate),
#                  'RMSprop':RMSprop(lr=learning_rate), 'Adadelta':Adadelta(lr=learning_rate),
#                  'Adagrad':Adagrad(lr=learning_rate), 'Adamax':Adamax(lr=learning_rate),
#                  'Nadam':Nadam(lr=learning_rate), 'Ftrl':Ftrl(lr=learning_rate)}
#     activationL = ['relu', 'sigmoid', 'softplus', 'softsign', 'tanh', 'selu',
#                    'elu', 'exponential', LeakyReLU,'relu']
#     neurons = round(neurons)
#     activation = activationL[round(activation)]
#     batch_size = round(batch_size)
#     epochs = round(epochs)
#     def nn_cl_fun():
#         opt = Adam(lr = learning_rate)
#         nn = Sequential()
#         nn.add(Dense(neurons, input_dim=X_train_scaled.shape[1], activation=activation))
#         nn.add(Dense(neurons, activation=activation))
#         nn.add(Dense(len(np.unique(labels)), activation='sigmoid'))
#         nn.compile(loss='sparse_categorical_crossentropy', optimizer=opt, metrics=['accuracy'])
#         return nn
#     es = EarlyStopping(monitor='accuracy', mode='max', verbose=0, patience=20)
#     nn = KerasClassifier(build_fn=nn_cl_fun, epochs=epochs, batch_size=batch_size,
#                          verbose=0)
#     kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)
#     score = cross_val_score(nn, X_train_scaled, y_train, scoring=score_acc, cv=kfold, fit_params={'callbacks':[es]}).mean()
#     return score

# # Set paramaters
# params_nn ={
#     'neurons': (10, 100),
#     'activation':(0, 9),
#     'optimizer':(0,7),
#     'learning_rate':(0.01, 1),
#     'batch_size':(200, 1000),
#     'epochs':(20, 100)
# }
# # Run Bayesian Optimization
# nn_bo = BayesianOptimization(nn_cl_bo, params_nn, random_state=111)
# nn_bo.maximize(init_points=25, n_iter=4)

# params_nn_ = nn_bo.max['params']
# activationL = ['relu', 'sigmoid', 'softplus', 'softsign', 'tanh', 'selu',
#                'elu', 'exponential', LeakyReLU,'relu']
# params_nn_['activation'] = activationL[round(params_nn_['activation'])]
# print(params_nn_)