import numpy as np
import random
import os
import os.path
import re
import argparse

from tensorflow.keras.models import Sequential
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras import optimizers
import tensorflow as tf
import env_modelue
from tensorflow.keras import layers, models
from tensorflow.keras import layers, Model
from tensorflow.keras.layers import Dense, Conv1D, MaxPooling1D, Dropout, Flatten, MultiHeadAttention

def CNN_only():
    model = Sequential()
    model.add(layers.Conv1D(input_shape=(66, 5), 
                            activation="relu", 
                            filters=100, 
                            kernel_size=10, 
                            strides=1, 
                            padding="same"))
    model.add(layers.MaxPooling1D(pool_size=2, strides=2))
    model.add(layers.Dropout(0.2))
    model.add(layers.Conv1D(activation="relu", 
                            filters=100, 
                            kernel_size=3, 
                            strides=1, 
                            padding="same"))
    model.add(layers.MaxPooling1D(pool_size=2, strides=2))
    model.add(layers.Dropout(0.2))
    model.add(layers.Flatten())
    model.add(layers.Dense(750, activation='relu'))
    model.add(layers.Dropout(0.2))
    model.add(layers.Dense(300, activation='relu'))
    model.add(layers.Dense(1, activation='sigmoid'))
    return model

def resnet_block(input_data, filters, conv_size):
    x = layers.Conv1D(filters, conv_size, activation='relu', padding='same')(input_data)
    x = layers.BatchNormalization()(x)
    x = layers.Conv1D(filters, conv_size, activation=None, padding='same')(x)
    x = layers.BatchNormalization()(x)
    x = layers.Add()([x, input_data])
    x = layers.Activation('relu')(x)
    x = layers.Dropout(0.2)(x) 
    return x

def resnet1d(input_shape=(66, 5), d_model=128, num_layers=6, num_heads=8, dff=512, dropout_rate=0.2):
    def res_block(x, filters, kernel_size=3, stride=1, dropout_rate=0.2):
        shortcut = x
        
        x = layers.Conv1D(filters, kernel_size, strides=stride, padding='same')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Activation('relu')(x)
        
        x = layers.Conv1D(filters, kernel_size, strides=1, padding='same')(x)
        x = layers.BatchNormalization()(x)
        
        if stride != 1:
            shortcut = layers.Conv1D(filters, 1, strides=stride, padding='same')(shortcut)
            shortcut = layers.BatchNormalization()(shortcut)
        
        x = layers.add([x, shortcut])
        x = layers.Activation('relu')(x)
        x = layers.Dropout(dropout_rate)(x)  # Dropout in each residual block
    
        return x
    
    inputs = tf.keras.Input(shape=input_shape)
    # Initial Conv Layer
    x = layers.Conv1D(128, 7, strides=2, padding='same')(inputs)
    x = layers.BatchNormalization()(x)
    x = layers.Activation('relu')(x)
    #x = res_block(x, 128, stride=2)
    x = res_block(x, 128, stride=1)
    #x = res_block(x, 128, stride=1)
    # Linear projection to d_model dimension

    
    x = layers.GlobalAveragePooling1D()(x)
    
    # Output
    outputs = layers.Dense(1, activation='sigmoid')(x)
    
    model = tf.keras.Model(inputs, outputs)
    return model
def DISMIR_deep():
    model = Sequential()
    model.add(layers.Conv1D(input_shape=(66, 5), 
                            activation="relu", 
                            filters=100, 
                            kernel_size=10, 
                            strides=1, 
                            padding="same"))
    model.add(layers.MaxPooling1D(pool_size=2, strides=2))
    model.add(layers.Dropout(0.2))
    model.add(layers.Bidirectional(layers.LSTM(33, return_sequences=True)))
    model.add(layers.Conv1D(input_shape=(33, 132), 
                            activation="relu", 
                            filters=100, 
                            kernel_size=3, 
                            strides=1, 
                            padding="same"))
    model.add(layers.MaxPooling1D(pool_size=2, strides=2))
    model.add(layers.Dropout(0.2))
    model.add(layers.Flatten())
    model.add(layers.Dense(750, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dropout(0.2))
    model.add(layers.Dense(300, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dense(1, activation='sigmoid', kernel_regularizer=None, bias_regularizer=None))
    return model

def lstm():
    model = Sequential()
    model.add(layers.Conv1D(input_shape=(66, 5), 
                            activation="relu", 
                            filters=100, 
                            kernel_size=10, 
                            strides=1, 
                            padding="same"))
    model.add(layers.MaxPooling1D(pool_size=2, strides=2))
    model.add(layers.Dropout(0.2))
    model.add(layers.Bidirectional(layers.LSTM(33, return_sequences=True)))
    model.add(layers.Dropout(0.2))
    model.add(layers.Flatten())
    model.add(layers.Dense(750, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dropout(0.2))
    model.add(layers.Dense(300, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dense(1, activation='sigmoid', kernel_regularizer=None, bias_regularizer=None))
    return model
def cnn_lstm():
    model = Sequential()
    model.add(layers.Bidirectional(layers.LSTM(33, return_sequences=True)))
    model.add(layers.Dropout(0.2))
    model.add(layers.Flatten())
    model.add(layers.Dense(750, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dropout(0.2))
    model.add(layers.Dense(300, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dense(1, activation='sigmoid', kernel_regularizer=None, bias_regularizer=None))
    return model

class PositionalEncoding(layers.Layer):
    def __init__(self, max_len, d_model):
        super(PositionalEncoding, self).__init__()
        self.max_len = max_len
        self.d_model = d_model
        
    def get_angles(self, pos, i, d_model):
        angle_rates = 1 / np.power(10000, (2 * (i//2)) / np.float32(d_model))
        return pos * angle_rates
    
    def call(self, inputs):
        seq_len = inputs.shape.as_list()[1]
        angle_rads = self.get_angles(np.arange(seq_len)[:, np.newaxis],
                                     np.arange(self.d_model)[np.newaxis, :],
                                     self.d_model)
        
        # apply sin to even indices in the array; 2i
        angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])
        
        # apply cos to odd indices in the array; 2i+1
        angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])
        
        pos_encoding = angle_rads[np.newaxis, ...]
        return inputs + tf.cast(pos_encoding, tf.float32)

def transformer_encoder(input_shape=(66, 5), d_model=128, num_layers=6, num_heads=8, dff=512, dropout_rate=0.1):
    inputs = tf.keras.Input(shape=input_shape)
    # Linear projection to d_model dimension
    x = layers.Dense(d_model)(inputs)
    # Add positional encoding
    x = PositionalEncoding(input_shape[0], d_model)(x)
    # Transformer layers
    for _ in range(num_layers):
        # Multi-head self-attention mechanism
        attn_output = layers.MultiHeadAttention(num_heads=num_heads, key_dim=d_model // num_heads)(x, x)
        attn_output = layers.Dropout(dropout_rate)(attn_output)
        # Residual connection and layer normalization
        out1 = layers.LayerNormalization(epsilon=1e-6)(x + attn_output)
        
        # Position-wise feed forward network
        ffn_output = layers.Dense(dff, activation='relu')(out1)
        ffn_output = layers.Dense(d_model)(ffn_output)
        ffn_output = layers.Dropout(dropout_rate)(ffn_output)
        # Residual connection and layer normalization
        x = layers.LayerNormalization(epsilon=1e-6)(out1 + ffn_output)
    
    
    x = layers.GlobalAveragePooling1D()(x)
    
    # Output
    outputs = layers.Dense(1, activation='sigmoid')(x)
    
    model = tf.keras.Model(inputs, outputs)
    return model
    
def res_transformer(input_shape=(66, 5), d_model=128, num_layers=6, num_heads=8, dff=512, dropout_rate=0.2):
    def res_block(x, filters, kernel_size=3, stride=1, dropout_rate=0.2):
        shortcut = x
        
        x = layers.Conv1D(filters, kernel_size, strides=stride, padding='same')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Activation('relu')(x)
        
        x = layers.Conv1D(filters, kernel_size, strides=1, padding='same')(x)
        x = layers.BatchNormalization()(x)
        
        if stride != 1:
            shortcut = layers.Conv1D(filters, 1, strides=stride, padding='same')(shortcut)
            shortcut = layers.BatchNormalization()(shortcut)
        
        x = layers.add([x, shortcut])
        x = layers.Activation('relu')(x)
        x = layers.Dropout(dropout_rate)(x)  # Dropout in each residual block
    
        return x
    
    inputs = tf.keras.Input(shape=input_shape)
    # Initial Conv Layer
    x = layers.Conv1D(128, 7, strides=2, padding='same')(inputs)
    x = layers.BatchNormalization()(x)
    x = layers.Activation('relu')(x)
    x = res_block(x, 128, stride=1)
    x = layers.Dense(d_model)(x)
    # Add positional encoding
    x = PositionalEncoding(input_shape[0], d_model)(x)
    
    # Transformer layers
    for _ in range(num_layers):
        # Multi-head self-attention mechanism
        attn_output = layers.MultiHeadAttention(num_heads=num_heads, key_dim=d_model // num_heads)(x, x)
        attn_output = layers.Dropout(dropout_rate)(attn_output)
        # Residual connection and layer normalization
        out1 = layers.LayerNormalization(epsilon=1e-6)(x + attn_output)
        
        # Position-wise feed forward network
        ffn_output = layers.Dense(dff, activation='relu')(out1)
        ffn_output = layers.Dense(d_model)(ffn_output)
        ffn_output = layers.Dropout(dropout_rate)(ffn_output)
        # Residual connection and layer normalization
        x = layers.LayerNormalization(epsilon=1e-6)(out1 + ffn_output)
    
    
    x = layers.GlobalAveragePooling1D()(x)
    
    # Output
    outputs = layers.Dense(1, activation='sigmoid')(x)
    
    model = tf.keras.Model(inputs, outputs)
    return model

def res_lstm(input_shape=(66, 5)):
    def res_block(x, filters, kernel_size=3, stride=1, dropout_rate=0.2):
        shortcut = x
        
        x = layers.Conv1D(filters, kernel_size, strides=stride, padding='same')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Activation('relu')(x)
        
        x = layers.Conv1D(filters, kernel_size, strides=1, padding='same')(x)
        x = layers.BatchNormalization()(x)
        
        if stride != 1:
            shortcut = layers.Conv1D(filters, 1, strides=stride, padding='same')(shortcut)
            shortcut = layers.BatchNormalization()(shortcut)
        
        x = layers.add([x, shortcut])
        x = layers.Activation('relu')(x)
        x = layers.Dropout(dropout_rate)(x)  # Dropout in each residual block
        
        return x
    
    inputs = tf.keras.Input(shape=input_shape)
    # Initial Conv Layer
    x = layers.Conv1D(128, 7, strides=2, padding='same')(inputs)
    x = layers.BatchNormalization()(x)
    x = layers.Activation('relu')(x)
    x = res_block(x, 128, stride=1)
    x = layers.Bidirectional(layers.LSTM(128))(x)
    x = layers.Dense(750, activation='relu')(x)
    x = layers.Dropout(0.2)(x)
    x = layers.Dense(300, activation='relu')(x)
    x = layers.Dense(1, activation='sigmoid')(x)
    #model.add(layers.Dense(1, activation='sigmoid', kernel_regularizer=None, bias_regularizer=None))
    # Output
    #outputs = layers.Dense(1, activation='sigmoid')(x)
    
    model = tf.keras.Model(inputs, x)
    return model

def cnn_transformer(input_shape=(66, 5), d_model=128, num_layers=6, num_heads=8, dff=512, dropout_rate=0.2):
    inputs = tf.keras.Input(shape=input_shape)
    x = layers.Conv1D(activation="relu", filters=100, kernel_size=10, strides=1, padding="same")(inputs)
    x = layers.MaxPooling1D(pool_size=2, strides=2)(x)
    # Add positional encoding
    x = layers.Dense(d_model)(x)
    x = PositionalEncoding(input_shape[0], d_model)(x)
    
    # Transformer layers
    for _ in range(num_layers):
        # Multi-head self-attention mechanism
        attn_output = layers.MultiHeadAttention(num_heads=num_heads, key_dim=d_model // num_heads)(x, x)
        attn_output = layers.Dropout(dropout_rate)(attn_output)
        # Residual connection and layer normalization
        out1 = layers.LayerNormalization(epsilon=1e-6)(x + attn_output)
        
        # Position-wise feed forward network
        ffn_output = layers.Dense(dff, activation='relu')(out1)
        ffn_output = layers.Dense(d_model)(ffn_output)
        ffn_output = layers.Dropout(dropout_rate)(ffn_output)
        # Residual connection and layer normalization
        x = layers.LayerNormalization(epsilon=1e-6)(out1 + ffn_output)
    
    
    x = layers.GlobalAveragePooling1D()(x)
    
    # Output
    outputs = layers.Dense(1, activation='sigmoid')(x)
    
    model = tf.keras.Model(inputs, outputs)
    return model


def cnn_lstm_transformer(input_shape=(66, 5), d_model=128, num_layers=6, num_heads=8, dff=512, dropout_rate=0.2):
    inputs = tf.keras.Input(shape=input_shape)
    x = layers.Conv1D(activation="relu", filters=100, kernel_size=10, strides=1, padding="same")(inputs)
    x = layers.MaxPooling1D(pool_size=2, strides=2)(x)
    # Add positional encoding
    x = layers.Dropout(0.2)(x)
    x = layers.Bidirectional(layers.LSTM(33, return_sequences=True))(x)
    x = layers.Dense(d_model)(x)
    x = PositionalEncoding(input_shape[0], d_model)(x)
    
    # Transformer layers
    for _ in range(num_layers):
        # Multi-head self-attention mechanism
        attn_output = layers.MultiHeadAttention(num_heads=num_heads, key_dim=d_model // num_heads)(x, x)
        attn_output = layers.Dropout(dropout_rate)(attn_output)
        # Residual connection and layer normalization
        out1 = layers.LayerNormalization(epsilon=1e-6)(x + attn_output)
        
        # Position-wise feed forward network
        ffn_output = layers.Dense(dff, activation='relu')(out1)
        ffn_output = layers.Dense(d_model)(ffn_output)
        ffn_output = layers.Dropout(dropout_rate)(ffn_output)
        # Residual connection and layer normalization
        x = layers.LayerNormalization(epsilon=1e-6)(out1 + ffn_output)
    
    
    x = layers.GlobalAveragePooling1D()(x)
    
    # Output
    outputs = layers.Dense(1, activation='sigmoid')(x)
    
    model = tf.keras.Model(inputs, outputs)
    return model


def cnn_lstm_cnn_transformer(input_shape=(66, 5), d_model=128, num_layers=6, num_heads=8, dff=512, dropout_rate=0.2):
    inputs = tf.keras.Input(shape=input_shape)
    x = layers.Conv1D(activation="relu", filters=100, kernel_size=10, strides=1, padding="same")(inputs)
    x = layers.MaxPooling1D(pool_size=2, strides=2)(x)
    # Add positional encoding
    x = layers.Dropout(0.2)(x)
    x = layers.Bidirectional(layers.LSTM(33, return_sequences=True))(x)
    x = layers.Conv1D(input_shape=(33, 132), activation="relu", filters=100, kernel_size=3, strides=1, padding="same")(x)
    x = layers.MaxPooling1D(pool_size=2, strides=2)(x)
    x = layers.Dense(d_model)(x)
    x = PositionalEncoding(input_shape[0], d_model)(x)
    
    # Transformer layers
    for _ in range(num_layers):
        # Multi-head self-attention mechanism
        attn_output = layers.MultiHeadAttention(num_heads=num_heads, key_dim=d_model // num_heads)(x, x)
        attn_output = layers.Dropout(dropout_rate)(attn_output)
        # Residual connection and layer normalization
        out1 = layers.LayerNormalization(epsilon=1e-6)(x + attn_output)
        
        # Position-wise feed forward network
        ffn_output = layers.Dense(dff, activation='relu')(out1)
        ffn_output = layers.Dense(d_model)(ffn_output)
        ffn_output = layers.Dropout(dropout_rate)(ffn_output)
        # Residual connection and layer normalization
        x = layers.LayerNormalization(epsilon=1e-6)(out1 + ffn_output)
    
    
    x = layers.GlobalAveragePooling1D()(x)
    x = layers.Dense(1, activation='sigmoid')(x)
    
    model = tf.keras.Model(inputs, x)
    return model

if __name__ == '__main__':
    x = tf.random.normal([32,66, 5])
    model = cnn_lstm()
    print(model(x).shape)
    model = lstm()
    print(model(x).shape)
    model = res_lstm()
    print(model(x).shape)
    model = cnn_transformer()
    print(model(x).shape)
    model = cnn_lstm_transformer()
    print(model(x).shape)
    model = cnn_lstm_cnn_transformer()
    print(model(x).shape)
    model = transformer_encoder()
    print(model(x).shape)