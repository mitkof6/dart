before_install() {
  cd /tmp
  PROJECTS='tinyxml2 assimp flann libccd fcl'
  for REPO in $PROJECTS
  do
    git clone git://github.com/dartsim/$REPO.git
    (cd $REPO; cmake .; make && sudo make install)
  done
  
  # Set github and bitbucket free of host checking questions
  echo -e "Host github.com\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config
  echo -e "Host bitbucket.org\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config
  # Install console_bridge
  git clone git://github.com/ros/console_bridge.git
  (cd console_bridge; cmake .; make && sudo make install)
  # Install urdfdom_headers
  hg clone https://bitbucket.org/osrf/urdfdom_headers
  (cd urdfdom_headers; cmake .; make && sudo make install)
  # Install console_bridge
  hg clone https://bitbucket.org/osrf/urdfdom
  (cd urdfdom; cmake .; make && sudo make install)
  # Install sdformat
  hg clone https://bitbucket.org/osrf/sdformat
  (cd sdformat; hg up sdf_1.4; cmake .; make && sudo make install)
}

APT='freeglut3 freeglut3-dev libglu1-mesa-dev libboost-all-dev cmake
cmake-curses-gui libeigen3-dev libxmu-dev libxi-dev libwxgtk2.8-dev
libgtest-dev libtinyxml-dev libgtest-dev'

sudo apt-get install $APT
(before_install)
