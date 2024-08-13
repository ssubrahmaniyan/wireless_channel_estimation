import json
import matplotlib.pyplot as plt
import numpy as np
with open("channels.txt", "r") as f:
    channel_values = [json.loads(line) for line in f]
channel_values = np.array(channel_values)

plt.plot(np.arange(0, channel_values[250:10000, 0].shape[0], 1), 10*np.log(abs(channel_values[250:10000, 0])))
#for i in range(1):
    #plt.plot(np.arange(0, channel_values[:500, i].shape[0], 1), ((channel_values[:500, i])))
#plt.plot(np.arange(0, channel_values[250:1000, 1].shape[0], 1), 10*np.log(abs(channel_values[250:1000, 1])), '--')
plt.show()