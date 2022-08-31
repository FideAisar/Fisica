---
tags: [Wiki]
---
 ## Enable and setup Bluetooth devices
---
 If you are using [pipewire](https://discovery.endeavouros.com/wp-admin/post.php?post=2519&action=edit) (default since Atlantis release):

`sudo pacman -S --needed bluez bluez-utils`

And graphical GUI tools to configure and manage Bluetooth:
- [blueberry](https://github.com/linuxmint/blueberry) (gtk) [recommended for XFCE]  
- [bluedevil](https://www.archlinux.org/packages/extra/x86_64/bluedevil/) (qt) [Integrate the Bluetooth technology within KDE workspace and applications]  
- [blueman](https://github.com/blueman-project/blueman) (gtk) [can be used in Desktop Environments independent]

DEs like KDE/plasma and GNOME ship their own Bluetooth GUI by default.

To use Bluetooth you will mainly need to enable it for use by default to start on every boot, or start it for a single session to use it only once:
- `sudo systemctl start bluetooth` –> This is only for single session use.
- `sudo systemctl enable bluetooth` –> Enables to start Bluetooth when booting the machine (permanent)

If you want to have Bluetooth enabled by default and start the service immediately:
`sudo systemctl enable --now bluetooth`