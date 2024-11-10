# 编译器和编译选项
CC = icc
CFLAGS = -Wall -Wextra -std=c11 -O0 -I$(INCLUDEDIR)

# 目标和目录
TARGET = build/w3test
SRCDIR = src
INCLUDEDIR = include

# 源文件和对象文件列表
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:.c=.o)

# 默认目标
all: $(TARGET)

# 链接目标文件生成可执行文件
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# 通用规则：编译每个.c文件为.o文件
$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# 清理规则
clean:
	rm -f $(SRCDIR)/*.o $(TARGET)

.PHONY: all clean
