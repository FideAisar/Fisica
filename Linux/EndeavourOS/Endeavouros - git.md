## Git setup
if you can't push to origin try:
```
git remote set-url origin $HTTPS/SSH link

sudo usermod -a -G $(stat -c '%G' .git) $USER
sudo chmod g+u .git -R
sudo chmod g+u .gitignore
su - $USER
```
