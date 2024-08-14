import json
import matplotlib.pyplot as plt
import numpy as np

with open("channels.txt", "r") as f:
    # Deserialize each line into complex numbers
    channel_values = [
        np.array([complex(re, im) for re, im in json.loads(line)])
        for line in f
    ]

# Convert the list of NumPy arrays to a 2D NumPy array
channel_values = np.array(channel_values)
channel_count = channel_values.shape[1]
for i in range(channel_count):
    plt.plot(np.arange(0, channel_values[250:10000, i].shape[0], 1), 10*np.log(abs(channel_values[250:10000, i])))
    #plt.plot(np.arange(0, channel_values[:500, i].shape[0], 1), ((channel_values[:500, i])))
#plt.plot(np.arange(0, channel_values[250:1000, 1].shape[0], 1), 10*np.log(abs(channel_values[250:1000, 1])), '--')
plt.show()