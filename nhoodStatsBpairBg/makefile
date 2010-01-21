include ${HOME}/kent/src/inc/common.mk
HG_INC += -I ${HOME}/kent/src/inc

L += -lm
MYLIBDIR = ${HOME}/kent/src/lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkweb.a

A = nhoodStatsBpairBg
O = nhoodStatsBpairBg.o

nhoodStatsBpairBg: ${O} ${MYLIBS}
	${CC} ${COPT} -o ${DESTDIR}${BINDIR}/${A}${EXE} $O ${MYLIBS} $L
	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

compile:: ${O}
	${CC} ${COPT} -o ${A}${EXE} ${O} ${MYLIBS} $L

clean::
	rm -f ${A} ${O}
