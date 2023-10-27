from tensorflow.keras import layers, models


def launch(input_layer, hidden_layers):
    output = input_layer
    for hidden_layer in hidden_layers:
        output = hidden_layer(output)
    return output


def model(length, kernel_size=10, filters=512, dense_ns=512):
    forward_input = layers.Input(shape=(length, 4))
    reverse_input = layers.Input(shape=(length, 4))
    hidden_layers = [
        layers.Conv1D(filters=filters, kernel_size=kernel_size),
        layers.LeakyReLU(alpha=0.1),
        layers.GlobalMaxPooling1D(),
        layers.Dropout(0.1),
    ]
    forward_output = launch(forward_input, hidden_layers)
    reverse_output = launch(reverse_input, hidden_layers)
    output = layers.Concatenate()([forward_output, reverse_output])
    output = layers.Dense(dense_ns, activation='relu')(output)
    output = layers.Dropout(0.1)(output)
    output = layers.Dense(2, activation='softmax')(output)
    model_ = models.Model(inputs=[forward_input, reverse_input], outputs=output)
    model_.compile(optimizer="adam", loss='binary_crossentropy', metrics='accuracy')
    return model_
