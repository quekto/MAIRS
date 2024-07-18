import matplotlib.pyplot as plt
import pandas as pd

csv_file_path = 'mairs/2024/0718_20K_CO_5min/0718OPow.csv'
df = pd.read_csv(csv_file_path)

x_name, y_name = df.columns[0], df.columns[1] # Get axis labels.

wavelength, spectra = df[x_name], df[y_name] # Get values

# Plot the data
plt.figure(figsize=(8, 6))
plt.plot(wavelength, spectra[::-1], marker='o', linestyle='-', color='b', label='Data')

# Adding title and labels
plt.title('X vs Y Plot')
plt.xlabel(x_name)
plt.ylabel(y_name)

# Adding grid
plt.grid(True)

# Adding legend
plt.legend()

# Show the plot
plt.show()
