# Компилятор
CXX = g++
CXXFLAGS = -Wall -O2 -std=c++17 -Wno-unused-variable

# Файлы
SRCS = f1.cpp f2.cpp riem_cuda.cpp
OBJS = $(SRCS:.cpp=.o)

# Исполняемый файл
EXEC = solver.exe

# Основной таргет: сборка всего проекта
all: $(EXEC)

# Компоновка исполняемого файла
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

# Компиляция каждого .cpp в .o
%.o: %.cpp func.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Очистка
clean:
	rm -f $(OBJS) $(EXEC)
