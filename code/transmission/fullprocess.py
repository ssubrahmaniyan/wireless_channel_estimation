import numpy as np
import librosa
import sounddevice as sd
from scipy.linalg import lstsq
import json
import soundfile as sf

audio_path = 'harvard.wav'
y, sr = librosa.load(audio_path)
print(f"Original audio shape: {y.shape}, Sample rate: {sr}")

print("Playing original audio...")
sd.play(y, sr)
sd.wait()


with open("channels.txt", "r") as f:
    channel_values = [json.loads(line) for line in f]

#Check for validity of channel value length
if len(channel_values) < len(y):
    print("Warning: Not enough channel values. Repeating the available values.")
    channel_values = (channel_values * (len(y) // len(channel_values) + 1))[:len(y)]

#Conversion to np array
channel_values = np.array(channel_values)

#AWGN level to be added
noise_level = 0.005

#Processing
for i in range(3):
    print(f"\nProcessing and playing version {i+1}:")
    
    scaled = y * channel_values[:len(y), i]
    print("Playing scaled audio...")
    sd.play(scaled, sr)
    sd.wait()
    
    noise = np.random.normal(0, noise_level, y.shape)
    noisy = scaled + noise
    print("Playing noisy audio...")
    sd.play(noisy, sr)
    sd.wait()
    
    #LS separation
    A = np.column_stack((scaled, noisy - scaled))
    coefficients, residuals, rank, s = lstsq(A, noisy)
    cleaned = coefficients[0] * scaled
    print("Playing cleaned audio...")
    sd.play(cleaned, sr)
    sd.wait()
    
    #File saving
    sf.write(f'scaled_audio_{i+1}.wav', scaled, sr)
    sf.write(f'noisy_audio_{i+1}.wav', noisy, sr)
    sf.write(f'cleaned_audio_{i+1}.wav', cleaned, sr)

print("\nAll audio files have been processed and saved.")
