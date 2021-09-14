ADD_EXECUTABLE(cs_demo1 cs_demo1.c)
TARGET_LINK_LIBRARIES(cs_demo1 cxsparse)

ADD_EXECUTABLE(cs_demo2 cs_demo2.c cs_demo.c cs_demo.h)
TARGET_LINK_LIBRARIES(cs_demo2 cxsparse)

ADD_EXECUTABLE(cs_demo3 cs_demo3.c cs_demo.c cs_demo.h)
TARGET_LINK_LIBRARIES(cs_demo3 cxsparse)

ADD_EXECUTABLE(cs_di_demo1 cs_di_demo1.c )
TARGET_LINK_LIBRARIES(cs_di_demo1 cxsparse)

ADD_EXECUTABLE(cs_di_demo2 cs_di_demo2.c cs_demo.c cs_demo.h)
TARGET_LINK_LIBRARIES(cs_di_demo2 cxsparse)

ADD_EXECUTABLE(cs_di_demo3 cs_di_demo3.c cs_demo.c cs_demo.h)
TARGET_LINK_LIBRARIES(cs_di_demo3 cxsparse)

