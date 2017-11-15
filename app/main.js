const {app, BrowserWindow} = require('electron')

let win;

//Function to create window and initiate default settings
function createWindow() {
    win = new BrowserWindow({
        width: 800,
        height: 800,
        backgroundColor: '#FFFFFF'
    });

    win.loadURL(`file://${__dirname}/dist/index.html`);

    win.on('closed', function () {
        win = null;
    });
}

//when app is ready call createWindow function
app.on('ready', createWindow);

//when all windows are closed, it closes app itself
app.on('window-all-closed', function(){

    //MacOs does not close app when all windows are not available
    if(process.platform !== 'darwin'){
        app.quit();
    }
});

app.on('activate', function () {

    if(win === null){
        createWindow();
    }
});