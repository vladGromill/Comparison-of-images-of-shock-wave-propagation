import numpy as np
import matplotlib.pyplot as plt
import re
import os
import sys
from pathlib import Path
from PIL import Image  # Pillow для быстрого сохранения изображений

def read_fld_file(filename):
    PATH = str(Path(os.getcwd(), filename))
    with open(PATH, 'r') as file:
        lines = file.readlines()

    # Читаем размеры сетки (I, J)
    i_match = re.search(r'I=(\d+)', lines[2])
    j_match = re.search(r'J=(\d+)', lines[2])
    if not i_match or not j_match:
        raise ValueError("Не удалось определить размеры сетки I и J")
    i_size = int(i_match.group(1))
    j_size = int(j_match.group(1))

    expected_size_nodes = i_size * j_size  # Количество узлов
    expected_size_cells = (i_size - 1) * (j_size - 1)  # Количество ячеек

    data = {}
    current_line = 5

    # Чтение X, Y (предполагаем, что они заданы в узлах сетки)
    for var in ['X', 'Y']:
        data[var] = []
        while len(data[var]) < expected_size_nodes:
            if current_line >= len(lines):
                raise IndexError(f"Недостаточно данных для {var}")
            numbers = list(map(float, lines[current_line].split()))
            data[var].extend(numbers)
            current_line += 1
        data[var] = np.array(data[var]).reshape((j_size, i_size))

    # Чтение переменных в центрах ячеек
    for var in ['D', 'U', 'V', 'P', 'Al']:
        data[var] = []
        while len(data[var]) < expected_size_cells:
            if current_line >= len(lines):
                raise IndexError(f"Недостаточно данных для {var}")
            numbers = list(map(float, lines[current_line].split()))
            data[var].extend(numbers)
            current_line += 1
        data[var] = np.array(data[var]).reshape((j_size - 1, i_size - 1))
    
    return data

def add_dark_noise(image, noise_level=0.06, dark_factor=0.45):
    """
    Добавляет затемняющий шум к изображению.
    - noise_level: уровень случайного шума
    - dark_factor: насколько затемнить изображение (0 - без изменений, 1 - почти чёрное)
    """
    noise = np.random.normal(-dark_factor, noise_level, image.shape)  # Смещённый в минус шум
    noisy_image = image + noise
    return np.clip(noisy_image, 0, 1)  # Ограничиваем диапазон значений от 0 до 1

def save_picture(data, path_to_save, filename, add_noise_flag=False):
    """
    Сохраняет изображение с использованием Pillow для ускорения.
    """
    image = data['D'].T
    if add_noise_flag:
        image = add_dark_noise(image)

    # Нормализация данных для изображения
    image_normalized = (image - image.min()) / (image.max() - image.min())  # Нормализация к [0, 1]
    image_normalized = (image_normalized * 255).astype(np.uint8)  # Преобразование в 8-битный формат

    # Создание и сохранение изображения с помощью Pillow
    img = Image.fromarray(image_normalized, mode='L')  # 'L' для grayscale
    img.save(str(Path(path_to_save, filename)))

def main():
    if len(sys.argv) != 2:
        print("Usage: python process_fld.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]  # Получаем имя файла из аргумента
    PATH_FOLDER = str(Path(os.getcwd()))
    PATH_TO_SAVE = str(Path(os.getcwd(), "results_pictures_dataset"))

    try:
        data = read_fld_file(filename)
        
        # Сохраняем обычное изображение
        save_picture(data, PATH_TO_SAVE, f"{filename.split('.')[0]}.png")
        
        # Сохраняем изображение с затемняющим шумом
        save_picture(data, PATH_TO_SAVE, f"{filename.split('.')[0]}_noisy.png", add_noise_flag=True)
        
        os.remove(str(Path(PATH_FOLDER, filename)))  # Удаляем исходный .fld файл
        #print(f"Processed and deleted {filename}")
    except Exception as e:
        print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    main()