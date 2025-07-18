# 编译器和基本选项
CC = icc
BASE_CFLAGS = -Wall -Wextra -std=c11 -I$(INCLUDEDIR)

# 目标和目录
TARGET = build/w3test
SRCDIR = src
INCLUDEDIR = include
SRCODIR = build

# 源文件和对象文件列表
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c, $(SRCODIR)/%.o, $(SRCS))

# 检查是否定义了DEBUG变量
ifdef DEBUG
    # 调试模式编译选项
    CFLAGS = $(BASE_CFLAGS) -g -O0 -DDEBUG
    BUILD_TYPE = debug
else
    # 发布模式编译选项
    CFLAGS = $(BASE_CFLAGS) -O2 -DNDEBUG
    BUILD_TYPE = release
endif

# 默认目标
all: $(TARGET)
	@echo "Build completed in $(BUILD_TYPE) mode"

# 调试模式目标
debug:
	@$(MAKE) DEBUG=1 all

# 发布模式目标
release:
	@$(MAKE) all

# 链接目标文件生成可执行文件
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET) in $(BUILD_TYPE) mode..."
	$(CC) $(CFLAGS) -o $@ $^

# 通用规则：编译每个.c文件为.o文件
$(SRCODIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(SRCODIR)
	@echo "Compiling $< ($(BUILD_TYPE) mode)..."
	$(CC) $(CFLAGS) -c $< -o $@

# 清理规则
clean:
	@echo "Cleaning build files..."
	rm -rf $(SRCODIR)/*.o $(TARGET)

# 深度清理（包括目录）
distclean: clean
	rm -rf $(SRCODIR)

# 显示编译信息
info:
	@echo "Compiler: $(CC)"
	@echo "Build type: $(BUILD_TYPE)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "Sources: $(SRCS)"
	@echo "Objects: $(OBJS)"
	@echo "Target: $(TARGET)"

# 重新构建
rebuild: clean all

# 帮助信息
help:
	@echo "Available targets:"
	@echo "  all      - Build in release mode (default)"
	@echo "  debug    - Build in debug mode with -DDEBUG"
	@echo "  release  - Build in release mode with -DNDEBUG"
	@echo "  clean    - Remove object files and executable"
	@echo "  distclean- Remove all build artifacts"
	@echo "  rebuild  - Clean and build"
	@echo "  info     - Show build configuration"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make debug    # Build with DEBUG enabled"
	@echo "  make release  # Build with optimizations"
	@echo "  make DEBUG=1  # Alternative way to enable debug"

.PHONY: all debug release clean distclean info rebuild help