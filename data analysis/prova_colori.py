colors = [
    "#00008B",  # Dark Blue
    "#0000CD",  # Medium Blue
    "#4169E1",  # Royal Blue
    "#1E90FF",  # Dodger Blue
    "#00BFFF",  # Deep Sky Blue
    "#87CEEB",  # Sky Blue
    "#87CEFA",  # Light Sky Blue
    "#EEE8AA",  # Pale Goldenrod
    "#F0E68C",  # Khaki
    "#FFD700"   # Gold
]

# If you want to plot them, you can use matplotlib as follows:

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 2))
for i, color in enumerate(colors):
    plt.fill_between([i, i+1], 0, 1, color=color)
plt.xlim(0, 10)
plt.axis('off')
plt.show()