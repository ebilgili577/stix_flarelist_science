from typing import List
from tensorflow import keras


def build_model(
    input_dim: int,
    hidden_dims: List[int],
    learning_rate: float = 1e-3
) -> keras.Model:
    """
    Build and compile the localization neural network.
    
    Architecture: input_dim -> hidden_dims[0] -> ... -> hidden_dims[-1] -> 2
    
    Args:
        input_dim: Number of input features
        hidden_dims: List of hidden layer sizes (e.g., [1024, 512, 256, 128, 64])
        learning_rate: Learning rate for Adam optimizer
        
    Returns:
        Compiled Keras model
    """
    layers = [keras.layers.Input(shape=(input_dim,))]
    
    # Add hidden layers
    for dim in hidden_dims:
        layers.append(keras.layers.Dense(dim, activation="relu"))
    
    # Add output layer
    layers.append(keras.layers.Dense(2))  # Output: (x, y) coordinates
    
    model = keras.Sequential(layers)
    
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=learning_rate),
        loss="mse",
        metrics=["mae"],
    )
    
    return model
