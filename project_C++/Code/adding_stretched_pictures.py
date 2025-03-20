import os
import random
from tqdm import tqdm
from PIL import Image
from pathlib import Path

folder_path = Path(os.getcwd()) / 'results_pictures_dataset'  # 'results_pictures_dataset'

# Получаем список всех пар картинок
images = [f for f in os.listdir(folder_path) if f.endswith('.png')]

# Выбираем пары картинок (чистая и шумная)
clean_images = [f for f in images if '_noisy' not in f]  # Чистые изображения
noisy_images = [f for f in images if '_noisy' in f]  # Шумные изображения

# Делаем выборку 30-40% пар изображений
selected_pairs = random.sample(list(zip(clean_images, noisy_images)), int(len(clean_images) * random.uniform(0.3, 0.4)))

# Функция для растяжения изображения
def stretch_image(image_path, stretch_factor=2):
    """
    Растягивает изображение по ширине с сохранением высоты.
    :param image_path: Путь к изображению.
    :param stretch_factor: Коэффициент растяжения (новая ширина = 128 * stretch_factor).
    :return: Растянутое изображение.
    """
    img = Image.open(image_path)
    width, height = img.size

    new_width = int(width * stretch_factor)
    new_height = height  # Высота остается прежней

    # Растягиваем изображение
    img_stretched = img.resize((new_width, new_height), Image.LANCZOS)
    return img_stretched

# Аугментация и сохранение новых изображений
for clean_img_name, noisy_img_name in tqdm(selected_pairs):
    clean_img_path = folder_path / clean_img_name
    noisy_img_path = folder_path / noisy_img_name

    # Выбираем случайный коэффициент растяжения
    stretch_factor = random.choice([2, 2.5, 3, 3.5])  # Пример коэффициентов

    # Растягиваем чистое и шумное изображения
    stretched_clean_img = stretch_image(clean_img_path, stretch_factor)
    stretched_noisy_img = stretch_image(noisy_img_path, stretch_factor)

    stretched_clean_img_name = clean_img_name.replace('.png', f'_stretched_{stretch_factor}x.png')
    stretched_noisy_img_name = noisy_img_name.replace('.png', f'_noisy_stretched_{stretch_factor}x.png')

    stretched_clean_img.save(folder_path / stretched_clean_img_name)
    stretched_noisy_img.save(folder_path / stretched_noisy_img_name)

    # print(f'Аугментировано и сохранено: {stretched_clean_img_name}, {stretched_noisy_img_name}')