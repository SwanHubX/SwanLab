from swanlab.integration.keras import SwanLabLogger
import tensorflow as tf
import swanlab

# Initialize SwanLab
swanlab.init(
    project="keras_mnist",
    experiment_name="mnist_example",
    description="Keras MNIST Example"
    )

# Load and preprocess MNIST data
(x_train, y_train), (x_test, y_test) = tf.keras.datasets.mnist.load_data()
x_train = x_train.reshape(-1, 28, 28, 1).astype('float32') / 255.0
x_test = x_test.reshape(-1, 28, 28, 1).astype('float32') / 255.0

# Build a simple CNN model
model = tf.keras.Sequential([
    tf.keras.layers.Conv2D(32, 3, activation='relu', input_shape=(28, 28, 1)),
    tf.keras.layers.MaxPooling2D(),
    tf.keras.layers.Conv2D(64, 3, activation='relu'),
    tf.keras.layers.MaxPooling2D(),
    tf.keras.layers.Flatten(),
    tf.keras.layers.Dense(64, activation='relu'),
    tf.keras.layers.Dense(10, activation='softmax')
])

# Compile the model
model.compile(
    optimizer='adam',
    loss='sparse_categorical_crossentropy',
    metrics=['accuracy']
)

# Train the model with SwanLabLogger
model.fit(
    x_train, 
    y_train,
    epochs=5,
    validation_data=(x_test, y_test),
    callbacks=[SwanLabLogger()]
)
