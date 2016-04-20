/**
 * Created by jon on 07/04/16.
 */

function setCookie(cname, cvalue, exdays) {
    var d = new Date();
    d.setTime(d.getTime() + (exdays*24*60*60*1000));
    var expires = "expires="+d.toUTCString();
    document.cookie = cname + "=" + cvalue + "; " + expires + "; path=/";
}

function getCookie(cname) {
    var name = cname + "=";
    var ca = document.cookie.split(';');
    for(var i=0; i<ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0)==' ') c = c.substring(1);
        if (c.indexOf(name) == 0) return c.substring(name.length, c.length);
    }
    return "";
}

function handleInfoPanelState(time) {
    var state = getCookie('infoPanelState');
    if (state == "true") {
        $(".left").hide(time);
        $("#min").hide(time);
        $("#max").show(time);
    }

    else {
        $(".left").show(time);
        $("#max").hide(time);
        $("#min").show(time);
    }
}

function handleInfoPanelMinimize() {
    setCookie("infoPanelState", "true", 50);
    handleInfoPanelState(500);
}

function handleInfoPanelMaximise() {
    setCookie("infoPanelState", "false", 50);
    handleInfoPanelState(500);
}

$(document).ready(function() {
    $("#min").click(handleInfoPanelMinimize);
    $("#max").click(handleInfoPanelMaximise);
    handleInfoPanelState(0);
    $(".divider").show();
});