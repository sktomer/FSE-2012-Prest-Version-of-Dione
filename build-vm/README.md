## Build your own Virtual Machine

### Here's how you can spin up a Virtual Machine for PREST:

This folder contains files and scripts used to build a virtual machine using Vagrant and Oracle Virtual Box.  
This VM will contain the "PREST" Tool.

A vagrant tutorial can be found [here](https://docs.vagrantup.com/v2/getting-started/index.html)  
Install [vagrant] (https://www.vagrantup.com/downloads.html) and [virtualbox](https://www.virtualbox.org/wiki/Downloads) (preferrably latest versions) on host machine.
The binaries and documentation for Oracle Virtual Box can be found [here](https://www.virtualbox.org/). 
Deploys Base Vagrant Box from the link for [Ubuntu 14.04 Desktop](https://vagrantcloud.com/box-cutter/boxes/ubuntu1404-desktop).

Follow this process to create the tool VM from scratch:  
* Download all files/folders under the "build-vm" folder.  
* Install Vagrant, if not already installed.  
* Install Oracle Virtual Box, if not already installed.  
* Open a command window.  
* Navigate to the location where the files were downloaded.  
* Type "vagrant up".  
* After built, VM GUI will be displayed.
There might be a prompt to restart. 
Please wait until the vagrant window returns to a prompt before rebooting.
The Virtual Machine also might reboot while it is being renamed.  

Note: The login credentials are:
* Username: vagrant
* Password: vagrant

Please do not interact with the VM until all of the activity in the command window has completed.  

