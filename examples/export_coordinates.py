#!/usr/bin/env python3
"""
Скрипт для експорту координат чорних дір з .dat файлу в txt формат
"""

import csv
import os
import sys

# Константи
year = 365.25636 * 24. * 3600.  # рік в секундах
msun = 1476.62504  # маса Сонця в метрах (G=c=1)

# Параметри системи OJ 287 (можна змінити)
M1 = 38e10  # Первинна СМЧД в сонячних масах
M2 = 20e8   # Вторинна ЧД в сонячних масах
m1 = M1 * msun
m2 = M2 * msun
M = m1 + m2  # Повна маса

# Вхідний файл (можна вказати як аргумент командного рядка)
if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    # За замовчуванням - останній файл
    input_file = "data/OJ287_precession_25-12-22-1314/OJ287_precession.dat"

# Вихідний файл
output_file = input_file.replace('.dat', '_coordinates.txt')

print(f"Reading: {input_file}")
print(f"Output: {output_file}")

# Читання даних
tt = []
xx1, yy1, zz1 = [], [], []
xx2, yy2, zz2 = [], [], []
rr = []
orbits_list = []

with open(input_file) as f:
    reader = csv.reader(f, delimiter=' ')
    for row in reader:
        tt.append(float(row[0]))
        xx1.append(float(row[1]))  # x1 велика ЧД
        yy1.append(float(row[2]))  # y1
        zz1.append(float(row[3]))  # z1
        xx2.append(float(row[4]))  # x2 мала ЧД
        yy2.append(float(row[5]))  # y2
        zz2.append(float(row[6]))  # z2
        rr.append(float(row[7]))   # відстань
        orbits_list.append(float(row[9]))  # номер орбіти

print(f"Read {len(tt)} data points")

# Запис у txt
with open(output_file, 'w') as f:
    f.write("# OJ 287 Binary Black Hole Coordinates\n")
    f.write(f"# M1 = {M1:.3e} solar masses (Primary SMBH)\n")
    f.write(f"# M2 = {M2:.3e} solar masses (Secondary BH)\n")
    f.write("# Coordinates normalized by total mass M = M1 + M2\n")
    f.write("# Time in years\n")
    f.write("#\n")
    f.write("# Columns:\n")
    f.write("# time[years]  x1/M  y1/M  z1/M  x2/M  y2/M  z2/M  separation/M  orbit\n")
    f.write("#" + "=" * 100 + "\n")
    
    # Прорідження для зменшення розміру (макс ~10000 точок)
    step = max(1, len(tt) // 10000)
    count = 0
    
    for i in range(0, len(tt), step):
        t_years = tt[i] / year
        x1 = xx1[i] / M
        y1 = yy1[i] / M
        z1 = zz1[i] / M
        x2 = xx2[i] / M
        y2 = yy2[i] / M
        z2 = zz2[i] / M
        r = rr[i] / M
        orb = orbits_list[i]
        
        f.write(f"{t_years:.6f}  {x1:.6f}  {y1:.6f}  {z1:.6f}  {x2:.6f}  {y2:.6f}  {z2:.6f}  {r:.6f}  {orb:.1f}\n")
        count += 1

print(f"Saved {count} data points to {output_file}")
print("Done!")

