import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';

// General imports
import { HomeComponent } from './layout/home/home.component';
import { FooterComponent } from './layout/footer/footer.component';
import { NavbarComponent } from './layout/navbar/navbar.component';
import { HeaderComponent } from './layout/header/header.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';
import { NgxPaginationModule } from 'ngx-pagination';

// GOG component imports
import { GogHolderComponent } from './gog-holder/gog-holder.component';
import { MapViewComponent } from './gog-holder/map-view/map-view.component';
import { ShowUsrMarkersComponent } from './gog-holder/show-usr-markers/show-usr-markers.component';
import { AddEditUsrMarkersComponent } from './gog-holder/add-edit-usr-markers/add-edit-usr-markers.component';

// GOG services imports
import { MarkersService } from './services/gog/markers/markers.service';

// API consumption imports
import { HttpClientModule } from '@angular/common/http';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';

@NgModule({
  declarations: [
    AppComponent,
    GogHolderComponent,
    HomeComponent,
    FooterComponent,
    NavbarComponent,
    HeaderComponent,
    NotFoundComponent,
    MapViewComponent,
    ShowUsrMarkersComponent,
    AddEditUsrMarkersComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    HttpClientModule,
    FormsModule,
    ReactiveFormsModule,
    NgxPaginationModule
  ],
  providers: [MarkersService],
  bootstrap: [AppComponent]
})
export class AppModule { }
