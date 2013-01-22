MODULE_NAME = yassi
all:
	gcc -shared -Wall -I/usr/include/python2.7/ -lpython2.7 -o $(MODULE_NAME).so $(MODULE_NAME).c

clean:
	rm -rf $(MODULE_NAME).c.*
	rm -rf $(MODULE_NAME).so
