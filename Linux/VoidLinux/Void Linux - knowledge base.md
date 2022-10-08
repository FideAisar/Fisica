---
title: Void Linux - knowledge base
author: Void Linux 
tags: [Wiki]
date: 05.09.2022
---

# Installation GUIDE
	sudo void-installer
## Install and update packages 
	xbps-install -Suv

## Firmware
	xbps-install linux-firmware-amd

## Locales
To enable a certain locale, un-comment or add the relevant lines in `/etc/default/libc-locales` and force-reconfigure the `glibc-locales` package.

Set  `LANG=xxxx`  in  `/etc/locale.conf`
## Sudo
To edit sudoers file

	sudo visudo

## Keymaps
Specifies which keymap to use for the Linux console. Available keymaps are listed in `/usr/share/kbd/keymaps`. For example:

	KEYMAP=it

## Font
Specifies which font to use for the Linux console. Available fonts are listed in `/usr/share/kbd/consolefonts`. For example:

	FONT=eurlatgr

To fix ugly font:

	ln -s /usr/share/fontconfig/conf.avail/70-no-bitmaps.conf /etc/fonts/conf.d/` `xbps-reconfigure -f fontconfig

## Sound setup (a)
	xbps-install -S alsa-utils alsa-plugins alsa-lib alsa-firmware
Add user to **audio** group → `sudo usermod -a -G audio <username>`

-   If audio doesn't work check:
    1.  `alsamixer` everything is not muted
    2.  Create _alsa-base.conf_  if not existed → `vim /etc/modprobe.d/alsa-base.conf` Add: **options snd-hda-intel index=1,0** in the file

## Bluetooth
	xbps-install bluez

## Install DE and base apps, here Kde, vim, firefox.
	xbps-install kde5 kde5-baseapps dbus vim firefox

## Enable services
	ln -sf /etc/sv/dbus /var/service

	ln -sf /etc/sv/tlp /var/service

	ln -sf /etc/sv/alsa /var/service

	ln -sf /etc/sv/bluetoothd /var/service

	ln -sf /etc/sv/elogind /var/service

	ln -sf /etc/sv/sddm /var/service

---
## Programs

Install:

| Program   | desc                                        |
| --------- | ------------------------------------------- |
| `dbus`    | *provide a system bus and/or a session bus* |
| `elogind` | *manage logins and power management*        |
| `vim`     | *text editor*                               |

## Video

Install video driver:

|        | Drivers | Wayland | Xorg |
| ------ | ------- | ------- | ---- |
| NVIDIA | `nouveau`        | `mesa-dri`        | `xf86-video-nouveau`     |
| AMD    |         |         |      |


## Graphical session
- Wayland: `Weston` / `Sway` / `Wayfire` / `Hikari` / `Cage` / `River`
- Xorg: `xorg`

Wayland compositors require `elogind` or `seatd` service.
Qt5-based applications require installing the `qt5-wayland` package and setting the environment variable `QT_QPA_PLATFORM=wayland-egl`. Some KDE specific applications also require installing the `kwayland` package. EFL-based applications require setting the environment variable `ELM_DISPLAY=wl`.

## Fonts*

To customize font display in your graphical session, you can use configurations provided in `/usr/share/fontconfig/conf.avail/`. To do so, create a symlink to the relevant `.conf` file in `/etc/fonts/conf.d/`, then use `xbps-reconfigure` to reconfigure the `fontconfig` package.

To disable bitmap fonts:
```bash
ln -s /usr/share/fontconfig/conf.avail/70-no-bitmaps.conf /etc/fonts/conf.d/ 
xbps-reconfigure -f fontconfig
```

## Icons*
By default the system uses Adwaita theme for icons. If you want to customize it:
Install `gtk+3` package, specify other themes in `/etc/gtk-3.0/settings`.

## Audio*
1. Install `pipewire` package and `libspa-bluetooth` and `bluez` (for bluetooth to work).
2. Enable `pipewire` and `bluetoothd` packages

```bash
xbps-install -S pipewire libspa-bluetooth bluez
ln -sf /etc/sv/pipewire /var/service
ln -sf /etc/sv/bluetoothd /var/service
```
## Services
```bash
ln -sf /etc/sv/dbus /var/service
```