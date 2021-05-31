ToolDAQPath=ToolDAQ

CPPFLAGS= -Wno-reorder -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -Wl,--no-as-needed

CC=g++ -std=c++1y -g -fPIC -shared $(CPPFLAGS)
CCC= g++ -std=c++1y -g -fPIC  $(CPPFLAGS)


ZMQLib= -L $(ToolDAQPath)/zeromq-4.0.7/lib -lzmq 
ZMQInclude= -I $(ToolDAQPath)/zeromq-4.0.7/include/ 

BoostLib= -L $(ToolDAQPath)/boost_1_66_0/install/lib -lboost_date_time -lboost_serialization  -lboost_iostreams -lboost_system -lboost_filesystem
BoostInclude= -I $(ToolDAQPath)/boost_1_66_0/install/include

RootInclude=  -I $(ToolDAQPath)/root-6.06.08/install/include
 
WCSimLib= -L ToolDAQ/WCSimLib -lWCSimRoot
WCSimInclude= -I ToolDAQ/WCSimLib/include

GenieIncludeDir := $(shell genie-config --topsrcdir)
GenieInclude= -I$(GenieIncludeDir)/Framework -I$(GenieIncludeDir)
GenieLibs= `genie-config --libs` -lxml2
PythiaLibs= -L ToolDAQ/Pythia6Support/v6_424/lib -lPythia6
Log4CppLibs= -L ToolDAQ/log4cpp/lib -llog4cpp
Log4CppInclude= -I ToolDAQ/log4cpp/include/log4cpp

RATEventLib= -L ToolDAQ/RATEventLib/lib -lRATEvent
RATEventInclude= -I ToolDAQ/RATEventLib/include

MrdTrackLib= -L ToolDAQ/MrdTrackLib/src -lFindMrdTracks
MrdTrackInclude= -I ToolDAQ/MrdTrackLib/include


RootLib=  -L $(ToolDAQPath)/root-6.06.08/install/lib `root-config --glibs` -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -m64 -lGui -lGenVector -lMinuit -lGeom -lEG -lEGPythia6 -lEve #-lGL -lGLEW -lGLU

RawViewerLib= -L UserTools/PlotWaveforms -lRawViewer


DataModelInclude = $(RootInclude)
DataModelLib = $(RootLib)


MyToolsInclude =  $(RootInclude) `python3.6-config --cflags` $(MrdTrackInclude) $(WCSimInclude) $(RATEventInclude) $(GenieInclude) $(Log4CppInclude)
MyToolsLib = -lcurl $(RootLib) `python3.6-config --libs` $(MrdTrackLib) $(WCSimLib) $(RATEventLib) $(RawViewerLib) $(GenieLibs) $(PythiaLibs) $(Log4CppLibs)


all: lib/libStore.so lib/libLogging.so lib/libDataModel.so include/Tool.h lib/libMyTools.so lib/libServiceDiscovery.so lib/libToolChain.so Analyse RemoteControl NodeDaemon

Analyse: src/main.cpp | lib/libMyTools.so lib/libStore.so lib/libLogging.so lib/libToolChain.so lib/libDataModel.so lib/libServiceDiscovery.so
	@echo -e "\n*************** Making " $@ "****************"
	g++ -std=c++1y -g -fPIC $(CPPFLAGS) src/main.cpp -o Analyse -I include -L lib -lStore -lMyTools -lToolChain -lDataModel -lLogging -lServiceDiscovery -lpthread $(DataModelInclude) $(DataModelLib) $(MyToolsInclude)  $(MyToolsLib) $(ZMQLib) $(ZMQInclude)  $(BoostLib) $(BoostInclude)


lib/libStore.so: $(ToolDAQPath)/ToolDAQFramework/src/Store/*
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libStore.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Store/*.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libStore.so lib/        
	#$(CC)  -I include $(ToolDAQPath)/ToolDAQFramework/src/Store/*.cpp -o lib/libStore.so $(BoostLib) $(BoostInclude)


include/Tool.h: $(ToolDAQPath)/ToolDAQFramework/src/Tool/Tool.h
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Tool/Tool.h include/
	cp UserTools/*.h include/
	cp UserTools/*/*.h include/
	cp DataModel/*.h include/


lib/libToolChain.so: $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/* | lib/libLogging.so lib/libStore.so lib/libServiceDiscovery.so lib/libLogging.so lib/libDataModel.so
	@echo -e "/n*************** Making " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/UserTools/Factory/*.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/*.h include/
	$(CC) $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/ToolChain.cpp -I include -lpthread -L lib -lStore -lDataModel -lServiceDiscovery -lLogging -o lib/libToolChain.so $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude) $(MyToolsInclude)  $(BoostLib) $(BoostInclude)


clean: 
	@echo -e "\n*************** Cleaning up ****************"
	rm -f include/*.h
	rm -f lib/*.so
	rm -f Analyse
	rm -f UserTools/*/*.o
	rm -f DataModel/*.o

lib/libDataModel.so: DataModel/* lib/libLogging.so lib/libStore.so $(patsubst DataModel/%.cpp, DataModel/%.o, $(wildcard DataModel/*.cpp))
	@echo -e "\n*************** Making " $@ "****************"
	cp -f DataModel/*.h include/
	$(CC) DataModel/*.o -I include -L lib -lStore  -lLogging  -o lib/libDataModel.so $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude)  $(BoostLib) $(BoostInclude)

lib/libMyTools.so: UserTools/*/* UserTools/* include/Tool.h lib/libLogging.so lib/libStore.so $(patsubst UserTools/%.cpp, UserTools/%.o, $(wildcard UserTools/*/*.cpp)) |lib/libDataModel.so lib/libToolChain.so 
	@echo -e "\n*************** Making " $@ "****************"
	cp -f UserTools/*/*.h include/
	cp -f UserTools/*.h include/
	#$(CC)  UserTools/Factory/Factory.cpp -I include -L lib -lStore -lDataModel -lLogging -o lib/libMyTools.so $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)
	$(CC) UserTools/*/*.o -I include -L lib -lStore -lDataModel -lLogging -o lib/libMyTools.so $(MyToolsInclude) $(DataModelInclude) $(MyToolsLib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)

RemoteControl:
	cd $(ToolDAQPath)/ToolDAQFramework/ && make RemoteControl
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/RemoteControl ./

NodeDaemon:
	cd $(ToolDAQPath)/ToolDAQFramework/ && make NodeDaemon
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/NodeDaemon ./

lib/libServiceDiscovery.so: $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/* | lib/libStore.so
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libServiceDiscovery.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/ServiceDiscovery.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libServiceDiscovery.so lib/
	#$(CC) -I include $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/ServiceDiscovery.cpp -o lib/libServiceDiscovery.so -L lib/ -lStore  $(ZMQInclude) $(ZMQLib) $(BoostLib) $(BoostInclude)

lib/libLogging.so: $(ToolDAQPath)/ToolDAQFramework/src/Logging/* | lib/libStore.so
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libLogging.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Logging/Logging.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libLogging.so lib/
	#$(CC) -I include $(ToolDAQPath)/ToolDAQFramework/src/Logging/Logging.cpp -o lib/libLogging.so -L lib/ -lStore $(ZMQInclude) $(ZMQLib) $(BoostLib) $(BoostInclude)

update:
	@echo -e "\n*************** Updating ****************"
	cd $(ToolDAQPath)/ToolDAQFramework; git pull
	cd $(ToolDAQPath)/zeromq-4.0.7; git pull
	cd $(ToolDAQPath)/MrdTrackLib; git checkout . ; git pull; make -f Makefile.FNAL
	cd $(ToolDAQPath)/WCSimLib; git checkout . ; git pull; make
	cd $(ToolDAQPath)/RATEventLib; git checkout . ; git pull; make
	git pull

UserTools/%.o: UserTools/%.cpp lib/libStore.so include/Tool.h lib/libLogging.so lib/libDataModel.so lib/libToolChain.so
	@echo -e "\n*************** Making " $@ "****************"
	cp $(shell dirname $<)/*.h include
	-$(CCC) -c -o $@ $< -I include -L lib -lStore -lDataModel -lLogging $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)

target: remove $(patsubst %.cpp, %.o, $(wildcard UserTools/$(TOOL)/*.cpp))

remove:
	echo "removing"
	-rm UserTools/$(TOOL)/*.o

DataModel/%.o: DataModel/%.cpp lib/libLogging.so lib/libStore.so
	@echo -e "\n*************** Making " $@ "****************"
	cp $(shell dirname $<)/*.h include
	-$(CCC) -c -o $@ $< -I include -L lib -lStore -lLogging  $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)
